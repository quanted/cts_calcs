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
        self.urlStruct = "/episuiteapi/rest/episuite/{}/estimated"  # new way (cgi server 1)
        # self.urlStruct = "/rest/episuite/{}/estimated"  # old way (local machine)
        self.methods = None
        self.melting_point = 0.0
        self.epi_props = ['melting_point', 'boiling_point', 'water_solubility', 'vapor_pressure', 'henrys_law_constant', 'log_kow', 'log_koc']
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
                'resultKey': 'waterSolMgLEstimated',
                'methods': {'wskownt': "WSKOWNT", 'waternt': "WATERNT"}
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
        return {'structure': "", 'melting_point': None}

    
    def makeDataRequest(self, structure, calc, prop, method=None):
        _post = self.getPostData(calc, prop)
        _post['structure'] = structure
        if self.melting_point != None:
            _post['melting_point'] = self.melting_point

        logging.info("getting url...")

        _url = self.baseUrl + self.getUrl(prop)

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
            if not key == 'nodes':
                _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict, 'method': None})

        try:
            _filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
            logging.info("EPI Filtered SMILES: {}".format(_filtered_smiles))
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            _response_dict.update({'data': "Cannot filter SMILES for EPI data"})
            return _response_dict

        try:
            if request_dict.get('prop') == 'water_sol' or request_dict.get('prop') == 'vapor_press':                
                self.melting_point = self.getMeltingPoint(_filtered_smiles, request_dict.get('sessionid'))
            else:
                self.melting_point = None
            _result_obj = self.makeDataRequest(_filtered_smiles, request_dict['calc'], request_dict['prop']) # make call for data!

            _response_dict.update({'data': _result_obj})

            # # NOTE: EPI now returns 2 values for water solubility
            if request_dict.get('prop') == 'water_sol':
                logging.warning("water sol property..")
                for data_obj in _result_obj.get('data', {}):
                    # replacing method name with ALL CAPS version:
                    logging.warning("water sol data obj: {}".format(data_obj))
                    _method_names = self.propMap['water_sol']['methods']
                    logging.warning("methods: {}".format(_method_names))
                    data_obj['method'] = _method_names.get(data_obj['method'])
        
            return _response_dict

        except Exception as err:
            logging.warning("Exception occurred getting {} data: {}".format(err, request_dict['calc']))
            _response_dict.update({'data': "cannot reach {} calculator".format(request_dict['calc'])})
            logging.info("##### session id: {}".format(request_dict.get('sessionid')))
            return _response_dict


    def getMeltingPoint(self, structure, sessionid):
        """
        Gets mass of structure from Measured, tries
        TEST if not available in Measured. Returns 0.0
        if neither have mp value.
        """
        melting_point_request = {
            'calc': "measured",  # should prob be measured
            # 'props': ['melting_point'],
            'prop': 'melting_point',
            'chemical': structure,
            'sessionid': sessionid
        }
        # todo: catch measured errors, then try epi melting point..
        # request = NotDjangoRequest(melting_point_request)
        # melting_point_response = measured_views.request_manager(request)
        melting_point = 0.0
        # measured_mp_response = MeasuredCalc().data_request_handler(melting_point_request)
        # Ideally need to make request to measured worker/queue:
        # tasks.measuredTask.apply_async(args=[melting_point_request], queue='measured')
        # measured_mp_resonse = tasks.measuredTask.apply(args[melting_point_request], queue='measured')

        # Probably best to call REST endpoint for melting point:
        # The rest call (cts_rest.py) will return one prop
        measured_mp_response = requests.post(
                                    os.environ.get('CTS_REST_SERVER') + '/cts/rest/measured/run', 
                                    data=json.dumps(melting_point_request), 
                                    allow_redirects=True,
                                    verify=False)


        logging.warning("MELTING POINT RESPONSE: {}".format(measured_mp_response.content))

        # # convert to python dict
        try:
            measured_results_object = json.loads(measured_mp_response.content)
            melting_point = measured_results_object['data']
            melting_point = float(melting_point['data'])
            # for data_obj in measured_results_object:
            #     if data_obj['prop'] == 'melting_point':
            #         melting_point = float(data_obj['data'])
            # melting_point = measured_mp_response['data']
        except Exception as e:
            logging.warning("Error in calculator_epi.py: {}".format(e))
            melting_point = 0.0

        if melting_point == 0.0 or not isinstance(melting_point, float):
            logging.warning("Trying to get MP from TEST..")
            try:
                melting_point_request['calc'] = 'test'
                # request = NotDjangoRequest(melting_point_request)
                # test_melting_point_response = test_views.request_manager(request)
                # test_mp_response = TestCalc().data_request_handler(melting_point_request)
                test_mp_response = requests.post(
                                    os.environ.get('CTS_REST_SERVER') + '/cts/rest/test/run', 
                                    data=json.dumps(melting_point_request), 
                                    allow_redirects=True,
                                    verify=False)
                logging.warning("TEST MP RESPONSE CONTENT: {}".format(test_mp_response))
                # melting_point = json.loads(test_melting_point_response.content)[0]['data']
                melting_point = json.loads(test_mp_response.content)['data']['data']
                logging.warning("TEST MP VALUE: {}".format(melting_point))
            except Exception as e:
                logging.warning("Error in calculator_epi.py: {}".format(e))
                melting_point = 0.0

            logging.warning("TEST MP TYPE: {}:".format(type(melting_point)))

            if not isinstance(melting_point, float):
                melting_point = 0.0
        # else:
        #     melting_point = melting_point_obj['data']

        logging.warning("MELTING POINT VALUE: {}".format(melting_point))

        return melting_point
