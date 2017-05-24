import requests
import json
import logging
import os
	
from calculator import Calculator
import smilesfilter

try:
	from cts_app.cts_calcs.calculator import Calculator
except ImportError as e:
	from cts_calcs.calculator import Calculator

headers = {'Content-Type': 'application/json'}


class MeasuredCalc(Calculator):
	"""
	Measured Calculator
	A single URL call returns data
	for all the properties: melting point,
	boiling point, vapor pressure, water solubility,
	log_kow, and henry's law constant
	"""

	def __init__(self):
		Calculator.__init__(self)

		self.postData = {"smiles" : ""}
		self.name = "measured"
		self.baseUrl = os.environ['CTS_EPI_SERVER']
		self.urlStruct = "/episuiteapi/rest/episuite/measured"  # new way
		# self.urlStruct = "/rest/episuite/measured"  # old way
		self.request_timeout = 10
		self.melting_point = 0.0

		# map workflow parameters to test
		self.propMap = {
			'melting_point': {
			   'result_key': 'melting_point'
			},
			'boiling_point': {
			   'result_key': 'boiling_point'
			},
			'water_sol': {
			   'result_key': 'water_solubility'
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
			}
		}

		self.result_structure = {
			'structure': '',
			'propertyname': '',
			'propertyvalue': None
		}

	def getPostData(self):
		return {"structure": ""}

	def makeDataRequest(self, structure):
		
		_post = self.getPostData()
		_post['structure'] = structure
		_url = self.baseUrl + self.urlStruct

		logging.info("Measured URL: {}".format(_url))

		# return self.request_logic(_url, _post)

		try:
			response = requests.post(_url, data=json.dumps(_post), headers=self.headers, timeout=self.request_timeout)
		except requests.exceptions.ConnectionError as ce:
			logging.warning("connection exception: {}".format(ce))
			raise 
		except requests.exceptions.Timeout as te:
			logging.warning("timeout exception: {}".format(te))
			raise 
		else:
			self.results = response
			return response

	def getPropertyValue(self, requested_property, response):
		"""
		Returns CTS data object for a requested
		property (cts format)
		"""
		# make sure property is in measured's format:
		if not requested_property in self.propMap.keys():
			# requested prop name doesn't match prop keys..
			raise KeyError(
				"requested property: {} for Measured data doesn't match Measured's property keys".format(
					requested_property))

		properties_dict = response['properties']
		data_obj = {
			'calc': "measured",
			'prop': requested_property
		}

		try:
			measured_requested_property = self.propMap[requested_property]['result_key']

			if measured_requested_property in properties_dict.keys():
				data_obj['data'] = properties_dict[measured_requested_property]['propertyvalue']
			else:
				data_obj['data'] = "property not available".format(measured_requested_property)

			return data_obj

		except Exception as err:
			logging.warning("Error at Measured Calc: {}".format(err))
			raise err

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
			return False

		# successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
		# json_obj = json.loads(response.content)

		# TODO: verify if blank data, finding the source of the empty water sol values...
		return True


	def data_request_handler(self, request_dict):

		_filtered_smiles = ''
		_response_dict = {}
		_measured_data = {}

		# fill any overlapping keys from request:
		for key in request_dict.keys():
			_response_dict[key] = request_dict.get(key)
		_response_dict.update({'request_post': request_dict, 'method': None})

		try:
			_filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
			logging.info("Measured Filtered SMILES: {}".format(_filtered_smiles))
		except Exception as err:
			logging.warning("Error filtering SMILES: {}".format(err))
			_response_dict.update({'data': "Cannot filter SMILES"})
			return _response_dict

		_retries = 3
		while _retries > 0:

			try:
				_response = self.makeDataRequest(_filtered_smiles) # make call for data!
				logging.info("Response from Measured: {}".format(_response))
				_measured_data = json.loads(_response.content)
				logging.info("Measured Data: {}".format(_measured_data))
			except Exception as e:
				logging.warning("Exception making request to Measured: {}".format(e))
				_response_dict.update({'data': "data not found"})
				# return _response_dict

			logging.info("Measured Data: {}".format(_measured_data))

			try:
				_data_obj = self.getPropertyValue(request_dict['prop'], _measured_data)
				logging.info("data object: {}".format(_data_obj))
				_response_dict.update(_data_obj)
				return _response_dict
			except Exception as err:
				logging.warning("Exception occurred getting Measured data: {}".format(err))
				_response_dict.update({'data': "cannot reach {} calculator".format(request_dict['calc'])})

			logging.info("Retrying Measured request..")
			_retries = _retries - 1