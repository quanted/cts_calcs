import requests
import json
import logging
import os
import math
from .calculator import Calculator
from .chemical_information import SMILESFilter



class OperaCalc(Calculator):
    """
	OPERA Suite Calculator
	"""
    def __init__(self):
        Calculator.__init__(self)
        self.method = None
        self.postData = {"smiles" : ""}
        self.name = "opera"
        self.baseUrl = os.environ['CTS_OPERA_SERVER']
        self.urlStruct = "/opera/rest/run"
        self.request_timeout = 300  # 3 min timeout for OPERAWS
        self.props = ['kow_no_ph', 'melting_point', 'boiling_point', 'henrys_law_con', 'vapor_press', 'water_sol', 'ion_con', 'kow_wph', 'log_bcf', 'koc']
        self.opera_props = ['LogP_pred', 'MP_pred', 'BP_pred', 'LogVP_pred', 'LogWS_pred', 'pKa_a_pred',
            'pKa_b_pred', 'LogD55_pred', 'LogD74_pred', 'LogBCF_pred', 'LogKoc_pred']
        self.propMap = {
            'kow_no_ph': {
                'result_key': "LogP_pred",  # is this correct?
                'methods': None
            },
            'melting_point': {
                'result_key': "MP_pred",
                'methods': None
            },
            'boiling_point': {
                'result_key': "BP_pred",
                'methods': None
            },
            'vapor_press': {
                'result_key': "LogVP_pred",
                'methods': None
            },
            'henrys_law_con': {
                'result_key': "LogHL_pred",
                'methods': None
            },
            'water_sol': {
                'result_key': "LogWS_pred",
                'methods': None
            },
            'ion_con': {
                'result_key': ["pKa_a_pred", "pKa_b_pred"],
                'methods': {'pKa_a_pred': "pKa", 'pKa_b_pred': "pKb"}
            },
            'kow_wph': {
                'result_key': ["LogD55_pred", "LogD74_pred"],
                'methods': {'LogD55_pred': "LogD55", 'LogD74_pred': "LogD74"}
            },
            'log_bcf': {
                'result_key': "LogBCF_pred",
                'methods': None
            },
            'koc': {
                'result_key': "LogKoc_pred",
                'methods': None
            }
        }
        self.meta_info = {
            'metaInfo': {
                'model': "opera",
                'collection': "qed",
                'modelVersion': "2.3",
                'description': "OPERA is a free and open-source/open-data suite of QSAR models providing predictions on physicochemical properties, environmental fate and toxcicity endpoints as well as additional information including applicability domain and accuracy assessment.",
                'status': '',
                'timestamp': self.gen_jid(),
                'url': "https://github.com/kmansouri/OPERA",
                'props': ['melting_point', 'boiling_point', 'water_sol', 'vapor_press', 'ion_con', 'henrys_law_con', 'kow_no_ph', 'koc', 'log_bcf', 'kow_wph']
            }
        }

    def match_chemical_with_node(self, chem_to_match, nodes_list):
        """
        Gets 'node' data from 'nodes' key for a given chemical.
        """
        if not nodes_list or len(nodes_list) < 1:
            return False
        chem_list = []
        chem_to_match = chem_to_match.replace("\n", "").replace("\r", "")

        for node in nodes_list:
            node['chemical'] = node['chemical'].replace("\n", "").replace("\r", "")

            if node['chemical'] == chem_to_match:
                return node
            elif node['smiles'] == chem_to_match:
                return node

        return False

    def convert_units_for_cts(self, prop, data_obj):
        """
        Converts certain OPERA properties to units used by CTS.
            + prop - property name to convert.
            + data_obj - response/data to be sent back to user.
        """
        if prop != 'ion_con' and math.isnan(float(data_obj['data'])):
            return "NaN"
        if prop in ['vapor_press', 'henrys_law_con']:
            # Converts from log:
            data_obj['data'] = 10**float(data_obj['data'])
        elif prop == 'water_sol':
            # Converts log(mol/L) --> mg/L:
            data_obj['mass'] = data_obj.get('mass')  # uses mass for WS conversion
            data_obj['data'] = self.convert_water_solubility(data_obj)
        elif prop == 'ion_con':
            # Sets 'data' to "none" if "NaN":
            data_obj = self.check_ion_con_for_nan(data_obj)

        return data_obj['data']

    def convert_water_solubility(self, ws_data_obj):
        """
        Converts water solubility from log(mol/L) => mg/L.
        """
        _ws_result = None
        ws_data_val = float(ws_data_obj['data'])
        if isinstance(ws_data_obj['mass'], float) or isinstance(ws_data_obj['mass'], int):
            # _ws_result = 1000 * float(ws_data_obj['mass']) * 10**-(ws_data_val)
            _ws_result = (10**ws_data_val) * ws_data_obj['mass'] * 1000.0
        else:
            # Requests mass from Calculator
            json_obj = self.getMass({'chemical': ws_data_obj['chemical']})
            mass = json_obj['data'][0]['mass']
            # _ws_result = 1000 * mass * 10**-(ws_data_val)  # from TESTWS WS conversion
            _ws_result = (10**ws_data_val) * mass * 1000.0
        return _ws_result
    
    def parse_results_for_cts(self, response_dict, opera_results):
        """
        Parses OPERA results for CTS API and CTS websockets.
        """
        requested_props = response_dict.get('props')
        if not requested_props:
            requested_props = [response_dict.get('prop')]

        # todo: add no 'data' exception handling
        opera_results = opera_results['data']
        chem_nodes = response_dict.get('nodes')
        result_index = 0
        curated_list = []
        for smiles_data_obj in opera_results:
            for prop in requested_props:
                prop_name = self.propMap[prop]['result_key']  # gets opera prop name
                # curated_dict = dict(response_dict)  # sends all key:vals for each prop result
                curated_dict = {}
                curated_dict['prop'] = prop
                curated_dict['data'] = ""
                curated_dict['chemical'] = response_dict['chemical'][result_index]
                curated_dict['node'] = self.match_chemical_with_node(curated_dict['chemical'], chem_nodes)
                curated_dict['calc'] = "opera"

                if prop == 'kow_wph' and isinstance(prop_name, list):
                    ph = response_dict.get('ph', self.default_ph)
                    if float(ph) == 5.5:
                        curated_dict['data'] = smiles_data_obj[prop_name[0]]  # expecting prop_name=[5.5, 7.4]
                    elif float(ph) == 7.4:
                        curated_dict['data'] = smiles_data_obj[prop_name[1]]
                    else:
                        curated_dict['data'] = "N/A"
                    curated_list.append(curated_dict)
                elif prop == 'ion_con' and isinstance(prop_name, list):
                    # Handles props with multiple results/methods:
                    _curated_results = self.parse_prop_with_multi_results(prop, prop_name, smiles_data_obj, response_dict, curated_dict)
                    self.convert_units_for_cts(prop, _curated_results)
                    curated_list.append(_curated_results)
                else:
                    curated_dict['data'] = smiles_data_obj[prop_name]
                    curated_dict['data'] = self.convert_units_for_cts(prop, curated_dict)
                    curated_list.append(curated_dict)

            result_index += 1
        curated_list = self.remove_nodes_key(curated_list)
        return curated_list

    def parse_prop_with_multi_results(self, prop, prop_name, smiles_data_obj, response_dict, curated_dict):
        """
        Further parses any property with methods into individual
        data objects.
        """
        for _prop in prop_name:
            prop_label = self.propMap[prop]['methods'][_prop]
            curated_dict['data'] += "{}: {}\n".format(prop_label, smiles_data_obj[_prop])
        return curated_dict

    def check_ion_con_for_nan(self, curated_dict):
        """
        Checks pka values and returns "none" if
        both pka and pkb are NaN.
        """
        pkas = curated_dict['data'].split("\n")  # Does OPERA ever return > 1 pka/pkbs?
        pka = pkas[0].split(":")[1].replace(" ", "")
        pkb = pkas[1].split(":")[1].replace(" ", "")

        pka, pkb = round(float(pka), 2), round(float(pkb), 2)

        # Packages ion_con data like sparc and chemaxon:
        ion_con_dict = {'pKa': [], 'pKb': []}
        if math.isnan(pka) and math.isnan(pkb):
            curated_dict['data'] = "none"
            return curated_dict
        if not math.isnan(pka):
            ion_con_dict['pKa'] = [pka]
        if not math.isnan(pkb):
            ion_con_dict['pKb'] = [pkb]
        curated_dict['data'] = ion_con_dict
        return curated_dict

    def remove_nodes_key(self, results_list):
        """
        Removes nodes key from results.
        """
        new_results_list = []
        for result in results_list:
            if 'nodes' in result:
                del result['nodes']
            new_results_list.append(result)
        return new_results_list

    def remove_opera_db_duplicates(self, db_results):
        """
        Makes sure there aren't data duplicates for a chemical's
        opera p-chem db results. This is inspired by vapor_pressure having
        duplicate documents in the pchem collection.
        TODO: Re-run DB script to remove vapor_press duplicates (warning: slow).
        """
        num_vp_docs = 0
        vp_indices = []
        prop_index = 0
        for result in db_results:
            if result.get('prop') == 'vapor_press':
                num_vp_docs += 1
                vp_indices.append(prop_index)
            prop_index += 1
        if num_vp_docs == 2:
            del db_results[vp_indices[0]]  # removes a vp duplicate entry
        return db_results

    def curate_logd(self, db_results, requested_dict, ph):
        if not 'kow_wph' in requested_dict.get('props'):
            return db_results
        new_results = []
        for result in db_results:
            if result.get('prop') == 'kow_wph' and result.get('ph') == float(ph):
                new_results.append(result)
            elif result.get('prop') == 'kow_wph' and float(ph) != 5.5 and float(ph) != 7.4:
                result['data'] = "N/A"
                new_results.append(result)
            elif result.get('prop') != 'kow_wph':
                new_results.append(result)
        return new_results

    def makeDataRequest(self, smiles):
        _post = {'smiles': smiles}
        _url = self.baseUrl + self.urlStruct
        return self.request_logic(_url, _post)
    
    def request_logic(self, url, post_data):
        """
        Handles retries and validation of responses
        """
        _valid_result = False  # for retry logic
        _retries = 0
        while not _valid_result and _retries < self.max_retries:
            # retry data request to chemaxon server until max retries or a valid result is returned
            try:
                response = requests.post(url, data=json.dumps(post_data), headers=self.headers, timeout=self.request_timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    self.results = json.loads(response.content)
                    return self.results
                _retries += 1
            except Exception as e:
                logging.warning("Exception in calculator_opera.py: {}".format(e))
                _retries += 1
            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))
        self.results = "calc server not found"
        return self.results

    def validate_response(self, response):
        """
        Validates sparc response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("opera server response status: {}".format(response.status_code))
            logging.warning("opera server response: {}".format(response.content))
            return False
        return True

    def data_request_handler(self, request_dict):
        """
        Makes requests to the OPERA Suite server
        """

        OPERA_URL = os.environ.get("CTS_OPERA_SERVER")

        if not isinstance(request_dict.get('chemical'), list):
            request_dict['chemical'] = [request_dict['chemical']]

        _response_dict = {}

        # fill any overlapping keys from request:
        for key in request_dict.keys():
            _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict, 'method': None})

        try:
            _result_obj = self.makeDataRequest(request_dict['chemical'])
            _result_obj = self.parse_results_for_cts(_response_dict, _result_obj)
            _response_dict['data'] = _result_obj
            _response_dict['valid'] = True
            return _response_dict
        except Exception as err:
            logging.warning("Exception occurred getting {} data: {}".format(request_dict['calc'], err))
            _response_dict.update({
                'data': "Cannot reach OPERA calculator",
                'valid': False
            })
            return _response_dict