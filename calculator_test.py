import requests
import json
import logging
import os
#try:
#    from cts_app.cts_calcs.calculator import Calculator
#except ImportError as e:
#    from cts_calcs.calculator import Calculator
from calculator import Calculator
from smilesfilter import parseSmilesByCalculator

headers = {'Content-Type': 'application/json'}


class TestCalc(Calculator):
    """
	TEST Suite Calculator
	"""

    def __init__(self):

        Calculator.__init__(self)

        self.postData = {"smiles" : ""}
        self.name = "test"
        self.baseUrl = os.environ['CTS_TEST_SERVER']
        self.urlStruct = "/api/TEST/{}/{}" # method, property

        # self.methods = ['hierarchical']
        self.methods = ['FDAMethod', 'HierarchicalMethod']
        # map workflow parameters to test
        self.propMap = {
            'melting_point': {
               'urlKey': 'MeltingPoint'
            },
            'boiling_point': {
               'urlKey': 'BoilingPoint'
            },
            'water_sol': {
               'urlKey': 'WaterSolubility'
            },
            'vapor_press': {
               'urlKey': 'VaporPressure'
            }
            # 'henrys_law_con': ,
            # 'kow_no_ph': 
        }


    def getPostData(self, calc, prop, method=None):
        return {"identifiers": {"SMILES": ""}}


    def makeDataRequest(self, structure, calc, prop, method=None):
        post = self.getPostData(calc, prop)
        post['identifiers']['SMILES'] = structure # set smiles
        test_prop = self.propMap[prop]['urlKey'] # prop name TEST understands
        url = self.baseUrl + self.urlStruct.format('FDAMethod', test_prop)
        try:
            response = requests.post(url, data=json.dumps(post), headers=headers, timeout=60)
        except requests.exceptions.ConnectionError as ce:
            logging.info("connection exception: {}".format(ce))
            # return None
            raise
        except requests.exceptions.Timeout as te:
            logging.info("timeout exception: {}".format(te))
            # return None
            raise

        self.results = response
        return response


    def convertWaterSolubility(self, _response_dict):
        """
        Converts water solubility from log(mol/L) => mg/L
        """
    	_ws_result = None
    	if isinstance(_response_dict['mass'], float) or isinstance(_response_dict['mass'], int):
    	        _ws_result = 1000 * float(_response_dict['mass']) * 10**-(_response_dict['test_datum'])
    	else:
    		# request mass from Calculator
    		json_obj = self.getMass({'chemical': _response_dict['chemical']})
    		mass = json_obj['data'][0]['mass']
    		_response_dict.update({'mass': mass})
            _ws_result = 1000 * float(_response_dict['mass'] * 10**-(_response_dict['test_datum']))
    	_response_dict.update({'data': _ws_result})
        return _response_dict


    def data_request_handler(self, request_dict):
        

        _filtered_smiles = ''
        _response_dict = {}

        # fill any overlapping keys from request:
        for key in request_dict.keys():
            _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict, 'method': None})


        # filter smiles before sending to TEST:
        # ++++++++++++++++++++++++ smiles filtering!!! ++++++++++++++++++++
        try:
            _filtered_smiles = parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            _response_dict.update({'data': "Cannot filter SMILES for TEST data"})
            return _response_dict
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        logging.info("TEST Filtered SMILES: {}".format(_filtered_smiles))

        try:
            logging.info("Calling TEST for {} data...".format(request_dict['prop']))

            _response = self.makeDataRequest(_filtered_smiles, request_dict['calc'], request_dict['prop'])
            _response_obj = json.loads(_response.content)
            #_response_dict.update({'data': _response_obj})

            logging.info("TEST response data for {}: {}".format(request_dict['prop'], _response))

            # sometimes TEST data successfully returns but with an error:
            if _response.status_code != 200:
                _response_dict['data'] = "TEST could not process chemical"
            else:
                _test_data = _response_obj['properties'][self.propMap[request_dict['prop']]['urlKey']]
		_response_dict['data'] = _test_data
		logging.warning("~~~ TEST DATA: {}".format(_test_data))
		logging.warning("~~~ MASS: {}".format(request_dict['mass']))
                if _test_data == -9999:
                    _response_dict['data'] = "N/A"
                elif request_dict['prop'] == 'water_sol':
                    _response_dict = self.convertWaterSolubility(_response_dict) # update response dict data
                #else:
                #    _response_dict['data'] = _test_data

            return _response_dict

        except Exception as err:
            logging.warning("Exception occurred getting TEST data: {}".format(err))
            _response_dict.update({'data': "timed out", 'request_post': request_dict})
            logging.info("##### session id: {}".format(request_dict.get('sessionid')))
            return _response_dict
