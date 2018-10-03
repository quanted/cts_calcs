import requests
import json
import logging
import os
#try:
#    from cts_app.cts_calcs.calculator import Calculator
#except ImportError as e:
#    from cts_calcs.calculator import Calculator
from .calculator import Calculator
from .chemical_information import SMILESFilter

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
				_ws_result = 1000 * float(_response_dict['mass']) * 10**-(_response_dict['data'])
		else:
			# request mass from Calculator
			json_obj = self.getMass({'chemical': _response_dict['chemical']})
			mass = json_obj['data'][0]['mass']
			_response_dict.update({'mass': mass})
			_ws_result = 1000 * float(_response_dict['mass'] * 10**-(_response_dict['data']))
		_response_dict.update({'data': _ws_result})
		return _response_dict


	def data_request_handler(self, request_dict):
		

		_filtered_smiles = ''
		_response_dict = {}

		# fill any overlapping keys from request:
		for key in request_dict.keys():
			if not key == 'nodes':
				_response_dict[key] = request_dict.get(key)
		_response_dict.update({'request_post': request_dict, 'method': None})


		# filter smiles before sending to TEST:
		# ++++++++++++++++++++++++ smiles filtering!!! ++++++++++++++++++++
		try:
			_filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
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
				# logging.warning("~~~ MASS: {}".format(request_dict['mass']))
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









class TestWSCalc(Calculator):
	"""
	TEST WS Calculator at https://comptox.epa.gov/dashboard/web-test/
	Additional documentation: Todd Martin's User's Guide
	"""

	def __init__(self):

		Calculator.__init__(self)

		self.postData = {"smiles" : ""}
		self.name = "test"

		# Example Request: https://comptox.epa.gov/dashboard/web-test/MP?smiles=CCCC&method=nn
		self.baseUrl = "https://comptox.epa.gov/dashboard/web-test/{}"  # input = property type (see propMap)

		self.methods = ['hc', 'nn', 'gc']  # general property methods
		self.method = None
		self.bcf_method = "sm"

		# self.methods_map = {
		# 	'Hierarchical clustering': 'hc',
		# 	'FDA': 'fda',
		# 	'Group contribution': 'gc',
		# 	'Nearest neighbor': 'nn'
		# 	# sm - single mode returning error still
		# }

		# map workflow parameters to test
		self.propMap = {
			'melting_point': {
			   'urlKey': 'MP'
			},
			'boiling_point': {
			   'urlKey': 'BP'
			},
			'water_sol': {
			   'urlKey': 'WS'
			},
			'vapor_press': {
			   'urlKey': 'VP'
			},
			'log_bcf': {
				'urlKey': 'BCF'
			}
			# 'henrys_law_con': ,
			# 'kow_no_ph': 
		}

		self.result_keys = ['id', 'smiles', 'expValMass', 'expValMolarLog', 'predValMass',
			'predValMolarLog', 'massUnits', 'molarLogUnits']

		self.cts_testws_data_key = 'predValMass'

		# TESTWS API responses map:
		self.response_map = {
			# NOTE: MP TESTWS endpoint is only returning '*ValMass', but with 'massUnits'="*C"
			'melting_point': {
				'data_type': 'predValMass'
			},
			# NOTE: BP TESTWS endpoint is only returning '*ValMass', but with 'massUnits'="*C"
			'boiling_point': {
				'data_type': 'predValMass'
			},
			'water_sol': {
				'data_type': 'predValMass'
			},
			'vapor_press': {
				'data_type': 'predValMass'
			},
			'log_bcf': {
				'data_type': 'predValMolarLog'
			}
		}



	def convertWaterSolubility(self, _response_dict):
		"""
		Converts water solubility from log(mol/L) => mg/L.
		Expecting water sol data from TESTWS to have the following keys:
		"expValMolarLog", "expValMass","predValMolarLog","predValMass","molarLogUnits","massUnits"
		"""
		# Requests mass from Jchem:
		json_obj = self.getMass({'chemical': _response_dict['chemical']})
		mass = json_obj['data'][0]['mass']
		_response_dict.update({'mass': mass})
		_ws_result = 1000 * float(_response_dict['mass']) * 10**-(float(_response_dict['data']))
		_response_dict.update({'data': _ws_result})
		return _response_dict



	def makeDataRequest(self, structure, calc, prop, method):
		test_prop = self.propMap[prop]['urlKey'] # prop name TEST understands
		# url = self.baseUrl + self.urlStruct.format('FDAMethod', test_prop)
		# _url = self.baseUrl + test_prop
		_url = self.baseUrl.format(test_prop)
		_payload = {'smiles': structure, 'method': method}
		try:
			# response = requests.post(url, data=json.dumps(post), headers=headers, timeout=10)
			response = requests.get(_url, params=_payload, timeout=10)
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



	def data_request_handler(self, request_dict):		

		_filtered_smiles = ''
		_response_dict = {}

		# fill any overlapping keys from request:
		for key in request_dict.keys():
			if not key == 'nodes':
				_response_dict[key] = request_dict.get(key)
		_response_dict.update({'request_post': request_dict})
		# _response_dict.update({'request_post': {'service': "pchemprops"}})  # TODO: get rid of 'request_post' and double data


		# filter smiles before sending to TEST:
		# ++++++++++++++++++++++++ smiles filtering!!! ++++++++++++++++++++
		try:
			_filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict.get('chemical'), self.name) # call smilesfilter
		except Exception as err:
			logging.warning("Error filtering SMILES: {}".format(err))
			_response_dict.update({'data': "Cannot filter SMILES for TEST WS data"})
			return _response_dict
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		# logging.info("TEST WS Filtered SMILES: {}".format(_filtered_smiles))
		# logging.info("Calling TEST WS for {} data...".format(request_dict['prop']))

		if request_dict.get('method') and request_dict['method'] in self.methods + [self.bcf_method]:
			# Uses method provided in request to get data from TESTWS, otherwise uses default
			self.method = request_dict.get('method')
			# Make sure method name is all caps (it's an acronym):
			_response_dict['method'] = _response_dict.get('method').upper()

		_response = self.makeDataRequest(_filtered_smiles, self.name, request_dict.get('prop'), self.method)

		if _response.status_code != 200:
			_response_dict.update({'data': "Cannot reach TESTWS"})
			return _response_dict

		_response_obj = json.loads(_response.content)
		_test_data = _response_obj['predictions'][0]  # list of predictions (getting first because only one chemical comes back for GET requests)

		if 'error' in _test_data:
			_response_dict.update({'data': "Cannot parse SMILES"})
			return _response_dict

		# Gets response key for property:
		data_type = self.response_map[request_dict['prop']]['data_type']

		# Sets response data to property's data key (based on desired units)
		if _test_data.get(data_type):
			_response_dict['data'] = _test_data[data_type]
	
		# Returns "N/A" for data if there isn't any TESTWS data found:
		if not 'data' in _response_dict or not _response_dict.get('data'):
			_response_dict['data'] = "N/A"
			return _response_dict

		# Reformats TESTWS VP result, e.g., "3.14*10^-15" -> "3.14e-15":
		if request_dict['prop'] == 'vapor_press':
			_response_dict['data'] = self.convert_testws_scinot(_response_dict['data'])

		return _response_dict



	def convert_testws_scinot(self, pchem_data):
		"""
		Converts TESTWS scientific notation format.
		Ex: "3.14*10^-15" -> "3.14e-15"
		"""
		try:
			split_val = pchem_data.split("*")  # splits up number for reformatting
			n = split_val[0]  # gets float portion
			p = split_val[1].split("^")[1]  # gets power portion
			new_num = "{}e{}".format(n, p)
			return new_num
		except Exception as e:
			logging.warning("Failed trying to reformat TESTWS VP.. Returning as-is..")
			logging.warning("Exception: {}".format(e))
			return pchem_data