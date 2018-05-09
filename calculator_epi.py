import requests
import json
import logging
import os
from .calculator import Calculator
from .calculator_measured import MeasuredCalc
from .calculator_test import TestCalc
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
        self.urlStruct = "/rest/episuite/estimated"  # newest way - local
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
                'result_key': 'log_koc'
            },
            'log_bcf': {
                'result_key': 'log_bcf',
                'methods': {'regression': "REG", 'Arnot-Gobas': "A-G"}
            },
            'log_baf': {
                'result_key': 'log_baf',
                'methods': {'Arnot-Gobas': "A-G"}
            }
        }


    def getPostData(self, calc, prop, method=None):
        return {'structure': "", 'melting_point': None}

    
    def makeDataRequest(self, structure, calc):
        # _post = self.getPostData(calc, prop)
        _post = {'structure': structure}
        if self.melting_point != None:
            _post['melting_point'] = self.melting_point


        # _url = self.baseUrl + self.getUrl(prop)
        _url = self.baseUrl + self.urlStruct

        logging.info("EPI URL: {}".format(_url))
        logging.info("EPI POST: {}".format(_post))

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
                    # break
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

        # successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
        # json_obj = json.loads(response.content)

        # TODO: verify if blank data, finding the source of the empty water sol values...
        return True


    def get_mp_from_results(self, results):
        for data_obj in results['data']:
                if data_obj.get('prop') == 'melting_point':
                    logging.info("Found MP in EPI results..")
                    return float(data_obj['data'])
                    
        return None


    # def request_manager(request):
    def data_request_handler(self, request_dict):
        """
        less_simple_proxy takes a request and
        makes the proper call to the TEST web services.
        it relies on the epi_calculator to do such.
        input: {"calc": [calculator], "prop": [property]}
        output: returns data from TEST server
        """

        EPI_URL = os.environ.get("CTS_EPI_SERVER")

        _filtered_smiles = ''
        _response_dict = {}

        logging.info("Request dict to epi: {}".format(request_dict))

        # fill any overlapping keys from request:
        for key in request_dict.keys():
            if not key == 'nodes':
                _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict, 'method': None})

        try:
            _filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
            logging.info("EPI Filtered SMILES: {}".format(_filtered_smiles))
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            _response_dict.update({
                'data': "Cannot filter SMILES",
                'valid': False
            })
            return _response_dict

        try:

            _get_mp = request_dict.get('prop') == 'water_sol' or request_dict.get('prop') == 'vapor_press'
            
            if _get_mp:
                self.melting_point = self.get_melting_point(_filtered_smiles, request_dict.get('sessionid'), 'epi')
            else:
                self.melting_point = None

            _result_obj = self.makeDataRequest(_filtered_smiles, request_dict['calc']) # make call for data!

            if _get_mp and not self.melting_point:
                # MP not found from measured or test, getting from results,
                # and requesting data again with set MP..
                logging.info("Trying to get MP from EPI request now..")
                self.melting_point = self.get_mp_from_results(_result_obj)
                logging.info("Making request to EPI with MP: {}".format(self.melting_point))
                _result_obj = self.makeDataRequest(_filtered_smiles, request_dict['calc'])  # Make request using MP

            _response_dict.update(_result_obj)
            _response_dict['valid'] = True
        
            return _response_dict

        except Exception as err:
            logging.warning("Exception occurred getting {} data: {}".format(err, request_dict['calc']))
            _response_dict.update({
                'data': "cannot reach {} calculator".format(request_dict['calc']),
                'valid': False
            })
            return _response_dict