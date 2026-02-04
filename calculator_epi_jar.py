import requests
import json
import logging
import os

from .calculator import Calculator
from .chemical_information import SMILESFilter
# from .calculator_rdkit import RdkitCalc
from .hydrolysis_jar import Hydrolysis



class EpiCalcJar(Calculator):
	"""
	EPI Suite Calculator
	"""
	def __init__(self):
		Calculator.__init__(self)
		self.hydrolysis = Hydrolysis()
		self.method = None
		self.postData = {"smiles" : ""}
		self.name = "epi"
		self.baseUrl = os.environ['CTS_EPI_SERVER']
		self.methods = None
		self.melting_point = None
		# self.epi_props = ['melting_point', 'boiling_point', 'water_solubility', 'vapor_pressure', 'henrys_law_constant', 'log_kow', 'koc', 'log_bcf', 'log_baf']
		self.epi_props = ['melting_point', 'boiling_point', 'water_sol', 'vapor_press', 'henrys_law_con', 'kow_no_ph', 'koc', 'log_bcf', 'log_baf']
		self.props = ['melting_point', 'boiling_point', 'water_sol', 'vapor_press', 'henrys_law_con', 'kow_no_ph', 'koc', 'log_bcf', 'log_baf']
		self.propMap = {
			'melting_point': {
			   'result_key': 'meltingPoint'
			},
			'boiling_point': {
			   'result_key': 'boilingPoint'
			},
			'water_sol': {
			   'result_key': ['waterSolubilityFromLogKow', 'waterSolubilityFromWaterNt'],
			   'methods': {'WSKOW': "WSKOW", 'WATERNT': "WATERNT"}  # Maintains similar pattern as calculator_epi
			},
			'vapor_press': {
			   'result_key': 'vaporPressure'
			},
			'henrys_law_con': {
				'result_key': 'henrysLawConstant'
			},
			'kow_no_ph': {
				'result_key': 'logKow'
			},
			'koc': {
				'result_key': 'logKoc',
				'methods': {'KOW': "KOW", "MCI": "MCI"}
			},
			'log_bcf': {
				'result_key': 'logBioconcentrationFactor',
				'methods': {'REG': "REG", 'A-G': "A-G"}  # Maintains similar pattern as calculator_epi
			},
			'log_baf': {
				'result_key': 'logBioaccumulationFactor',
				'methods': {'A-G': "A-G"}  # Maintains similar pattern as calculator_epi
			},
			'qsar': {
				'result_key': 'qsar',
			}
		}

	
	def makeDataRequest(self, url, structure, calc=None):
		_post = {"smiles": structure}
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
				response = requests.get(url, params=post_data)
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
		logging.warning("get_mp_from_results results: {}".format(results))
		for data_obj in results['data']:
				if data_obj.get('prop') == 'melting_point':
					logging.info("Found MP in EPI results..")
					return float(data_obj['data'])
		return None


	def parse_api_results(self, results, _filtered_smiles):
		"""
		Parses pchem results from EPI API into format used
		by the CTS backend.
		"""

		parsed_data = {"data": []}

		# Extract data and convert keys
		for cts_prop, api_info in self.propMap.items():

			api_key = api_info["result_key"]
			methods = api_info.get("methods")

			if api_key == "qsar":
				continue

			data_obj = {
				"chemical": _filtered_smiles,
				"calc": "epi",
				"prop": cts_prop,
				"method": None,
				"data": None
			}

			if cts_prop == "water_sol":
				# Multiple props for some cts methods (e.g., water_sol)
				method_vals = list(methods.values())
				i = 0
				for api_prop in api_key:
					estimated_value = results[api_prop].get("estimatedValue", {}).get("value", None)
					new_item = dict(data_obj)
					new_item["method"] = method_vals[i]
					new_item["data"] = str(estimated_value)
					parsed_data["data"].append(new_item)
					i += 1
			elif cts_prop == "log_bcf":
				estimated_value =  results.get("bioconcentration", {}).get("logBioconcentrationFactor", None)
				data1 = dict(data_obj)
				data1["method"] = methods["REG"]
				data1["data"] = str(estimated_value)
				parsed_data["data"].append(data1)
				
				estimated_value =  results.get("bioconcentration", {}).get("arnotGobasBcfBafEstimates", {})[0].get("logBioconcentrationFactor", None)
				data2 = dict(data_obj)
				data2["method"] = methods["A-G"]
				data2["data"] = str(estimated_value)
				parsed_data["data"].append(data2)
				
			elif cts_prop == "log_baf":
				estimated_value =  results.get("bioconcentration", {}).get("logBioaccumulationFactor", None)
				new_item = dict(data_obj)
				new_item["method"] = methods["A-G"]
				new_item["data"] = str(estimated_value)
				parsed_data["data"].append(new_item)
			elif cts_prop == "koc":
				log_koc_data = results.get("logKoc", {}).get("estimatedValue", {}).get("model", {}).get("models")
				log_koc_vals = {model["name"]: model["correctedLogKoc"] for model in log_koc_data}

				data1 = dict(data_obj)
				data1["method"] = methods["MCI"]
				data1["data"] = log_koc_vals["MCI"]
				parsed_data["data"].append(data1)

				data2 = dict(data_obj)
				data2["method"] = methods["KOW"]
				data2["data"] = log_koc_vals["Kow"]
				parsed_data["data"].append(data2)
			else:
				estimated_value = results[api_key].get("estimatedValue", {}).get("value", None)
				logging.warning("estimated_value: {}".format(estimated_value))
				new_item = dict(data_obj)
				if "methods" in api_info:
					new_item["method"] = list(methods.values())[0]
				new_item["data"] = str(estimated_value)
				parsed_data["data"].append(new_item)

		return parsed_data


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

			_result_obj = self.hydrolysis.make_qsar_request(request_dict)

			# TODO: Account for not valid result_obj

			_response_dict['data'] = _result_obj
			_response_dict['valid'] = True

			return _response_dict

		try:

			_get_mp = request_dict.get('prop') == 'water_sol' or request_dict.get('prop') == 'vapor_press'
			
			if _get_mp:
				self.melting_point = self.get_melting_point(_filtered_smiles, 
										request_dict.get('sessionid'), self)
			else:
				self.melting_point = None

			logging.warning("Making data request\nURL: {}\nSMILES: {}\n".format(self.baseUrl, _filtered_smiles))

			_result_obj = self.makeDataRequest(self.baseUrl, _filtered_smiles, request_dict['calc']) # make call for data!

			logging.warning("makeDataRequest result: {}".format(_result_obj))

			_result_obj = self.parse_api_results(_result_obj, _filtered_smiles)

			logging.warning("parsed result: {}".format(_result_obj))


			################################################
			# # TODO: Update MP request for new EPI API.
			################################################
			# if _get_mp and not self.melting_point:
			#     # MP not found from measured or test, getting from results,
			#     # and requesting data again with set MP..
			#     self.melting_point = self.get_mp_from_results(_result_obj)
			#     _result_obj = self.makeDataRequest(self.baseUrl, _filtered_smiles, request_dict['calc'])  # Make request using MP




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