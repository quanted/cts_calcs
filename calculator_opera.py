import requests
import json
import logging
import os
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
        self.request_timeout = 180
        self.props = ['kow_no_ph', 'melting_point', 'boiling_point', 'vapor_press', 'water_sol', 'ion_con', 'kow_wph', 'log_bcf', 'koc']
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
                # 'result_key': ["LogD55_pred", "LogD74_pred"],  # is this correct?
                'result_key': "LogD74_pred",  # is data for 5.5 and 7.4 pH values?
                'methods': None
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

    # def parse_results_for_cts(self, opera_results, requested_props):
    def parse_results_for_cts(self, response_dict, opera_results):
        """
        Parses OPERA results for CTS API and CTS websockets.
        """
        requested_props = response_dict.get('props')
        if not requested_props:
            requested_props = [response_dict.get('prop')]

        # todo: add no 'data' exception handling
        opera_results = opera_results['data']

        curated_list = []
        for smiles_data_obj in opera_results:
            for prop in requested_props:
                prop_name = self.propMap[prop]['result_key']  # gets opera prop name
                if isinstance(prop_name, list):
                    # Handles props with multiple results/methods:
                    _curated_results = self.parse_prop_with_multi_results(prop, prop_name, smiles_data_obj, response_dict)
                    curated_list.append(_curated_results)
                else:
                    # curated_dict = {}
                    curated_dict = dict(response_dict)  # sends all key:vals for each prop result
                    curated_dict['prop'] = prop
                    curated_dict['data'] = smiles_data_obj[prop_name]
                    curated_dict['calc'] = "opera"
                    curated_list.append(curated_dict)
        return curated_list

    def parse_prop_with_multi_results(self, prop, prop_name, smiles_data_obj, response_dict):
        """
        Further parses any property with methods into individual
        data objects.
        """
        # curated_dict = {}
        curated_dict = dict(response_dict)
        curated_dict['prop'] = prop
        curated_dict['calc'] = "opera"
        curated_dict['data'] = ""
        for _prop in prop_name:
            prop_label = self.propMap[prop]['methods'][_prop]
            prop_data = smiles_data_obj[_prop]
            curated_dict['data'] += "{}: {}\n".format(prop_label, prop_data)
        return curated_dict

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
                logging.warning("Exception in calculator_epi.py: {}".format(e))
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

        _filtered_smiles = ''
        _response_dict = {}

        # fill any overlapping keys from request:
        for key in request_dict.keys():
            if not key == 'nodes':
                _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict, 'method': None})

        try:
            _filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            _response_dict.update({
                'data': "Cannot filter SMILES",
                'valid': False
            })
            return _response_dict

        try:
            _result_obj = self.makeDataRequest(_filtered_smiles)

            # requested_props = request_dict.get('props')
            # if not requested_props:
            #     requested_props = [request_dict.get('prop')]
            # _result_obj = self.parse_results_for_cts(_result_obj['data'], requested_props)

            _result_obj = self.parse_results_for_cts(_response_dict, _result_obj)


            # _response_dict.update(_result_obj)
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

    # def convert_property(self, )