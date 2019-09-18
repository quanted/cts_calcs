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
		self.baseUrl = os.environ.get('CTS_TEST_SERVER')
		if not self.baseUrl:
			self.baseUrl = "https://comptox.epa.gov/dashboard/web-test"

		# TESTWS with bulk prediction endpoint (09/2019):
		self.bulk_predict_endpoint = "http://webtest2.sciencedataexperts.com/bulk/predict"

		# self.methods = ['hc', 'nn', 'gc']  # general property methods
		self.methods = ['hc', 'nn', 'gc', 'sm']  # general property methods
		self.method = None
		# self.bcf_method = "sm"

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

		# POST for testwsV2 requesting all props and methods in one request:
		self.bulk_predict_post = {
			"format": "smiles",
			"endpoints": ["mp", "bp", "ws", "vp", "bcf"],
			"methods": ["hc","gc", "nn", "sm"],
			"reportTypes": [],
			"query": ""  # smiles goes here
		}

		self.pchem_response = {
			'chemical': "",
			'calc': "testws",
			'prop': "",
			'method': None,
			'data': None
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



	def make_bulk_request(self, structure):
		"""
		Makes request to TESTWS to get all props and methods using
		the bulk predict endpoint from v2.
		"""
		url = self.bulk_predict_endpoint
		post = dict(self.bulk_predict_post)
		post['query'] = structure
		try:
			response = requests.post(url, data=json.dumps(post), headers=headers, timeout=10)
		except requests.exceptions.ConnectionError as ce:
			logging.info("connection exception: {}".format(ce))
			raise Exception("Cannot connect")
		except requests.exceptions.Timeout as te:
			logging.info("timeout exception: {}".format(te))
			raise Exception("Timeout")
		except Exception as e:
			raise Exception("Request error")
		return response



	def makeDataRequest(self, structure, calc, prop, method):
		test_prop = self.propMap[prop]['urlKey'] # prop name TEST understands
		_url = self.baseUrl + "/{}".format(test_prop) 
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
		"""
		Request handler for TESTWS bulk predict endpoint.
		"""
		_filtered_smiles = ''
		_response_dict = {}

		# fill any overlapping keys from request:
		for key in request_dict.keys():
			if not key == 'nodes':
				_response_dict[key] = request_dict.get(key)
		_response_dict.update({'request_post': request_dict})


		# filter smiles before sending to TESTWS:
		# ++++++++++++++++++++++++ smiles filtering!!! ++++++++++++++++++++
		try:
			_filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict.get('chemical'), self.name) # call smilesfilter
		except Exception as err:
			logging.warning("Error filtering SMILES: {}".format(err))
			_response_dict.update({'data': "Cannot filter SMILES for TESTWS"})
			return _response_dict
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		try:
			response = self.make_bulk_request(_filtered_smiles)
		except Exception as e:
			logging.warning("Exception in calculator_test: {}".format(e))
			return self.build_error_response(_response_dict, "Cannot reach TESTWS")

		if response.status_code != 200:
			return self.build_error_response(_response_dict, "TESTWS server error")

		testws_response_obj = json.loads(response.content)

		try:
			_response_dict['data'] = self.build_response_object(_filtered_smiles, testws_response_obj)
			_response_dict['status'] = True
			return _response_dict
		except Exception as e:
			logging.warning("Exception in calculator_test: {}".format(e))
			return self.build_error_response(_response_dict, "TESTWS parse error")



	def build_error_response(self, response_obj, message):
		"""
		Builds error response object for each testws prop.
		"""
		error_response = {'status': False, 'data': []}
		for prop_obj in self.propMap:
			error_obj = dict(self.pchem_response)
			error_obj['chemical'] = response_obj['chemical']
			error_obj['prop'] = prop_obj['urlKey']
			error_obj['data'] = message
			error_response['data'].append(error_obj)
		return error_response



	def build_response_object(self, chemical, response_obj):
		"""
		Loops TESTWS bulk predict response and builds response object
		for CTS.
		"""
		if not 'predictions' in response_obj:
			raise Exception("Error parsing results")
		if not isinstance(response_obj.get('predictions'), list):
			raise Exception("Error parsing results")
		if 'error' in response_obj['predictions'][0]:
			raise Exception("Error requesting data")

		test_results = response_obj['predictions']  # single-chemical list of results
		cts_results = []  # list of curated test results

		for data_obj in test_results:
			response_dict = dict(self.pchem_response)
			response_dict['chemical'] = chemical
			response_dict['method'] = data_obj['method'].upper()
			response_dict['prop'] = self.get_cts_prop_name(data_obj['endpoint'])
			response_dict['data'] = self.get_testws_data(response_dict, data_obj)
			cts_results.append(response_dict)

		return cts_results



	def get_testws_data(self, response_obj, test_results_obj):
		"""
		Picks out data key and performs any neccessary unit conversions.
		"""
		# Gets response key for property:
		data_type = self.response_map[response_obj['prop']]['data_type']
		# Sets response data to property's data key (based on desired units)
		if test_results_obj.get(data_type):
			return test_results_obj[data_type]
		# Returns "N/A" for data if there isn't any TESTWS data found:
		if not 'data' in response_obj or not response_obj.get('data'):
			return "N/A"
		# Reformats TESTWS VP result, e.g., "3.14*10^-15" -> "3.14e-15":
		if response_obj['prop'] == 'vapor_press':
			return self.convert_testws_scinot(response_obj['data'])



	def get_cts_prop_name(self, testws_prop_name):
		"""
		Returns cts prop name from testws prop name.
		"""
		for cts_name, val in self.propMap.items():
			testws_name = val['urlKey']
			if testws_name == testws_prop_name:
				return cts_name
		


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