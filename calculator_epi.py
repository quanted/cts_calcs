import requests
import json
import logging
import os

from calculator import Calculator
import smilesfilter

# try:
#     # from cts_app.cts_calcs.calculator import Calculator
#     from cts_app.cts_calcs.smilesfilter import max_weight
#     # from cts_app.cts_calcs import smilesfilter
# except ImportError as e:
#     # from cts_calcs.calculator import Calculator
#     from cts_calcs.smilesfilter import max_weight
#     # from cts_calcs import smilesfilter

headers = {'Content-Type': 'application/json'}


class EpiCalc(Calculator):
    """
	EPI Suite Calculator
	"""
    def __init__(self):
        Calculator.__init__(self)
        self.postData = {"smiles" : ""}
        self.name = "epi"
        self.baseUrl = os.environ['CTS_EPI_SERVER']
        # self.urlStruct = "/episuiteapi/rest/episuite/{}/estimated"  # new way (cgi server 1)
        self.urlStruct = "/rest/episuite/{}/estimated"  # old way (local machine)
        self.methods = None
        self.props = ['melting_point', 'boiling_point', 'water_sol', 'vapor_press', 'henrys_law_con', 'kow_no_ph', 'koc']
        self.propMap = {
            'melting_point': {
                'urlKey': 'meltingPoint',
                'propKey': 'Melting Pt (deg C)(estimated)',
                'resultKey': 'meltingPtDegCEstimated'
            },
            'boiling_point': {
                'urlKey': 'boilingPoint',
                'propKey': '',
                'resultKey': 'boilingPtDegCEstimated'
            },
            'water_sol': {
                'urlKey': 'waterSolubility',
                'propKey': '',
                'resultKey': 'waterSolMgLEstimated'
            },
            'vapor_press': {
                'urlKey': 'vaporPressure',
                'propKey': '',
                'resultKey': 'vaporPressMmHgEstimated'
            },
            'henrys_law_con': {
                'urlKey': 'henrysLawConstant',
                'propKey': '',
                'resultKey': 'henryLcBondAtmM3Mole'
            },
            'kow_no_ph': {
                'urlKey': 'logKow',
                'propKey': '',
                'resultKey': 'logKowEstimate'
            },
            'koc': {
                'urlKey': 'soilAdsorptionCoefficientKoc',
                'propKey': '',
                'resultKey': 'soilAdsorptionCoefKoc'
            }
        }


    def getPostData(self, calc, prop, method=None):
        return {'structure': ""}

    
    def makeDataRequest(self, structure, calc, prop, method=None):
        _post = self.getPostData(calc, prop)
        _post['structure'] = structure

        logging.info("getting url...")

        _url = self.baseUrl + self.getUrl(prop)

        logging.info("EPI URL: {}".format(_url))

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
                    break
                _retries += 1
            except Exception as e:
                logging.warning("Exception in calculator_sparc.py: {}".format(e))
                _retries += 1

            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))
        return self.results


    def validate_response(self, response):
        """
        Validates sparc response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("epi server response status: {}".format(response.status_code))
            return False

        # successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
        # json_obj = json.loads(response.content)

        # TODO: verify if blank data, finding the source of the empty water sol values...
        return True





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

        # fill any overlapping keys from request:
        for key in request_dict.keys():
            _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict})

        try:
            _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
            logging.info("EPI Filtered SMILES: {}".format(_filtered_smiles))
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            _response_dict.update({'data': "Cannot filter SMILES for EPI data"})
            return _response_dict

        try:
            _result_obj = self.makeDataRequest(_filtered_smiles, request_dict['calc'], request_dict['prop']) # make call for data!

            if 'propertyvalue' in _result_obj:
                _response_dict.update({'data': _result_obj['propertyvalue']})
            else:
                _response_dict.update({'data': _result_obj})

            return _response_dict

        except Exception as err:
            logging.warning("Exception occurred getting {} data: {}".format(err, request_dict['calc']))
            _response_dict.update({'data': "cannot reach {} calculator".format(request_dict['calc'])})
            logging.info("##### session id: {}".format(request_dict['sessionid']))
            return _response_dict