import requests
import json
import logging
import os
from .calculator import Calculator
from .chemical_information import SMILESFilter



class EpiCalc(Calculator):
    """
    EPI Suite Calculator
    """
    def __init__(self):
        Calculator.__init__(self)
        self.method = None
        self.postData = {"smiles" : ""}
        self.name = "epi"
        self.baseUrl = os.environ['CTS_EPI_SERVER']
        # self.urlStruct = "/episuiteapi/rest/episuite/estimated"  # newest way - server
        # self.urlStruct = "/rest/episuite/estimated"  # newest way - local
        self.methods = None
        self.melting_point = None
        self.epi_props = ['melting_point', 'boiling_point', 'water_solubility', 'vapor_pressure', 'henrys_law_constant', 'log_kow', 'log_koc', 'log_bcf', 'log_baf']
        self.props = ['melting_point', 'boiling_point', 'water_sol', 'vapor_press', 'henrys_law_con', 'kow_no_ph', 'koc', 'log_bcf', 'log_baf']
        self.propMap = {
            'melting_point': {
               'result_key': 'melting_point'
            },
            'boiling_point': {
               'result_key': 'boiling_point'
            },
            'water_sol': {
               'result_key': 'water_solubility',
               'methods': {'WSKOW': "WSKOW", 'WATERNT': "WATERNT"}
            },
            'vapor_press': {
               'result_key': 'vapor_pressure'
            },
            'henrys_law_con': {
                'result_key': 'henrys_law_constant'
            },
            'kow_no_ph': {
                'result_key': 'log_kow'
            },
            'koc': {
                'result_key': 'log_koc',
                'methods': {'MCI': "MCI", 'Kow': "KOW"}
            },
            'log_bcf': {
                'result_key': 'log_bcf',
                'methods': {'regression': "REG", 'Arnot-Gobas': "A-G"}
            },
            'log_baf': {
                'result_key': 'log_baf',
                'methods': {'Arnot-Gobas': "A-G"}
            },
            'qsar': {
                'result_key': 'qsar',
            }
        }
        self.qsar_request_map = {
            'halogenated aliphatics: elimination': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (no adjacent x)': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (vicinal x)': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (geminal x)': 'hydrolysis/alkylhalide',
            'epoxide hydrolysis': 'hydrolysis/epoxide',
            'organophosphorus ester hydrolysis 1': 'hydrolysis/phosphate',  # Phosphate or Thiophosphate
            'organophosphorus ester hydrolysis 2': 'hydrolysis/phosphate',  # Phosphate or Thiophosphate
            'carboxylic acid ester hydrolysis': 'hydrolysis/ester',
            'anhydride hydrolysis': 'hydrolysis/anhydride',
            'carbamate hydrolysis': 'hydrolysis/carbamate'
        }
        self.cleaved_list = [
            'organophosphorus ester hydrolysis 1',
            'organophosphorus ester hydrolysis 2',
            'carboxylic acid ester hydrolysis',
            'anhydride hydrolysis',
            'carbamate hydrolysis'
        ]
        self.op_esters = [
            'organophosphorus ester hydrolysis 1',
            'organophosphorus ester hydrolysis 2'
        ]


    def getPostData(self, calc, prop, method=None):
        return {'structure': "", 'melting_point': None}

    
    def makeDataRequest(self, url, structure, calc=None):
        _post = {'structure': structure}
        if self.melting_point != None:
            _post['melting_point'] = self.melting_point
        return self.request_logic(url, _post)

    
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
            logging.warning("epi server response status: {}".format(response.status_code))
            logging.warning("epi server response: {}".format(response.content))
            return False
        return True



    def get_mp_from_results(self, results):
        for data_obj in results['data']:
                if data_obj.get('prop') == 'melting_point':
                    logging.info("Found MP in EPI results..")
                    return float(data_obj['data'])
        return None


    def round_half_life(self, value):
        upper_bound = 1e3
        lower_bound = 1e-1

        if type(value) != float:
            value = float(value)

        if abs(value) > upper_bound or abs(value) < lower_bound:
            return "{:.2e}".format(value)
        else:
            return round(value, 2)


    def make_qsar_request(self, request_dict):
        """
        Makes requests to epi suite for half-lives.
        """
        structure = request_dict.get('filtered_smiles')
        route = request_dict.get('route').lower()
        route_url = None
        num_sites = None
        unique_schemes_count = int(request_dict.get("uniqueSchemesCount"))

        if not route in list(self.qsar_request_map.keys()):
            raise Exception("Route not found.")

        route_url = self.qsar_request_map[route]
        url = self.baseUrl.replace('estimated', '') + route_url

        if route in self.cleaved_list and route in self.op_esters:
            logging.info("Route in cleaved list and an OP Ester.")
            num_sites = request_dict.get("chemical").count("P")  # gets count of "P" from parent smiles
        elif route in self.cleaved_list and not route in self.op_esters:
            logging.info("Route in cleaved list but not OP Ester.")
            num_sites = int(request_dict.get("productCount")) / 2
        else:
            logging.info("Route not in cleavd list.")
            num_sites = int(request_dict.get("productCount"))

        logging.info("Incoming request_dict for QSAR request: {}".format(request_dict))
        logging.info("Number of sites: {}".format(num_sites))
        logging.info("Request to EPI for half life:\nURL:{}\nStructure:{}".format(url, structure))

        # TODO: Account for OP Ester route where num_sites > 1 (cases C and D).

        response = requests.post(url, data=json.dumps({'structure': structure}), headers=self.headers)

        if response.status_code != 200:
            return {
                "error": "Error getting QSAR data from EPI Suite",
                "valid": False
            }

        try:
            response_obj = json.loads(response.content)
            # if not response_obj.get("data") or len(response_obj.get("data")) != 1:
            if not response_obj.get("data") or len(response_obj.get("data")) < 1:
                raise Exception("Half-life response does not have excepted 'data' key or size is unexpected.\nResponse: {}".format(response_obj))
            
            # NOTE: If OP Ester 1/2 and num_sites <= 1, pick specific half-life value from set.
            # Get Kb if OP Ester 1, and Ka/n for OP Ester 2

            logging.info("Half-life response data: {}".format(response_obj["data"]))

            halfLifeValue = None
            for data_obj in response_obj["data"]:
                if num_sites <= 1:
                    logging.info("Number of sites <= 1")
                    if route == self.op_esters[0] and data_obj["prop"] == "Kb":
                        halfLifeValue = self.round_half_life(data_obj["data"])
                        break
                    elif route == self.op_esters[1] and (data_obj["prop"] == "Ka" or data_obj["prop"] == "Kn"):
                        halfLifeValue = self.round_half_life(data_obj["data"])
                        break
                else:
                    logging.info("Number of sites > 1")
                    # Case C or D
                    if unique_schemes_count > 1:
                        logging.info("unique_schemes_count > 1")
                        # Case C
                        pass  # don't see case C example in excel file
                    else:
                        logging.info("unique_schemes_count <= 1")
                        # Case D
                        if "halogenated aliphatics" in route.lower():
                            # return N/A
                            logging.info("Case D halogenated aliphatics")
                            halfLifeValue = "N/A"
                            break
                        elif "epoxide" in route.lower() and (data_obj["prop"] == "Ka" or data_obj["prop"] == "Kn"):
                            # return Ka/n
                            logging.info("Case D epoxide, returning Ka/n")
                            halfLifeValue = self.round_half_life(data_obj["data"])
                            break
                        
                        if route in self.cleaved_list:
                            if "anhydride" in route.lower():
                                logging.info("Case D cleaved anhydride")
                                # return kb1 kb2 factored in (on epi wrapper side)
                                halfLifeValue = self.round_half_life(data_obj["data"])
                                break
                            elif data_obj["prop"] == "Kb":
                                logging.info("Case D not anhydride, returning Kb")
                                # return Kb
                                halfLifeValue = self.round_half_life(data_obj["data"])
                                break

                        elif data_obj["prop"] == "Kb":
                            logging.info("Case D not cleaved, returning Kb")
                            # return Kb
                            halfLifeValue = self.round_half_life(data_obj["data"])
                            break

            if not halfLifeValue:
                halfLifeValue = self.round_half_life(response_obj["data"][0]["data"])

            response_obj["data"][0]["data"] = halfLifeValue

            return {
                "data": response_obj,
                "valid": True
            }
        except Exception as e:
            logging.warning("Error parsing QSAR response: {}".format(e))
            return {
                "error": "Error getting QSAR data from EPI Suite",
                "valid": False
            }


    def data_request_handler(self, request_dict):
        """
        Makes requests to the EPI Suite server
        """
        
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

        # Handle QSAR request or continue to the usual p-chem stuff
        if request_dict.get('prop') == 'qsar':
            
            request_dict['filtered_smiles'] = _filtered_smiles

            _result_obj = self.make_qsar_request(request_dict)

            # TODO: Account for not valid result_obj

            _response_dict.update(_result_obj['data'])
            _response_dict['valid'] = True

            return _response_dict

        try:

            _get_mp = request_dict.get('prop') == 'water_sol' or request_dict.get('prop') == 'vapor_press'
            
            if _get_mp:
                self.melting_point = self.get_melting_point(_filtered_smiles, 
                                        request_dict.get('sessionid'), self)
            else:
                self.melting_point = None

            _result_obj = self.makeDataRequest(self.baseUrl, _filtered_smiles, request_dict['calc']) # make call for data!

            if _get_mp and not self.melting_point:
                # MP not found from measured or test, getting from results,
                # and requesting data again with set MP..
                self.melting_point = self.get_mp_from_results(_result_obj)
                _result_obj = self.makeDataRequest(self.baseUrl, _filtered_smiles, request_dict['calc'])  # Make request using MP

            _response_dict.update(_result_obj)
            _response_dict['valid'] = True
        
            return _response_dict

        except Exception as err:
            logging.warning("Exception occurred getting {} data: {}".format(err, request_dict['calc']))
            _response_dict.update({
                'data': "Cannot reach EPI calculator",
                'valid': False
            })
            return _response_dict
