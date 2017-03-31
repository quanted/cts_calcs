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
		# self.urlStruct = "/episuiteapi/rest/episuite/measured"  # new way
		self.urlStruct = "/rest/episuite/measured"  # old way

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
		
		post = self.getPostData()
		post['structure'] = structure
		url = self.baseUrl + self.urlStruct

		logging.info("Measured URL: {}".format(url))

		try:
			response = requests.post(url, data=json.dumps(post), headers=headers, timeout=30)
		except requests.exceptions.ConnectionError as ce:
			logging.info("connection exception: {}".format(ce))
			raise 
		except requests.exceptions.Timeout as te:
			logging.info("timeout exception: {}".format(te))
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


	def data_request_handler(self, request_dict):

		_filtered_smiles = ''
		_response_dict = {}
		_measured_data = {}

		# fill any overlapping keys from request:
		for key in request_dict.keys():
			_response_dict[key] = request_dict.get(key)
		_response_dict.update({'request_post': request_dict})

		try:
			_filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
			logging.info("Measured Filtered SMILES: {}".format(_filtered_smiles))
		except Exception as err:
			logging.warning("Error filtering SMILES: {}".format(err))
			_response_dict.update({'data': "Cannot filter SMILES"})
			for prop in request_dict['props']:
				_response_dict.update({'prop': prop})
				self.redis_conn.publish(request_dict['sessionid'], json.dumps(_response_dict))
			return


		if request_dict['run_type'] == 'rest':
			request_dict['props'] = [request_dict['prop']]


		try:
			_response = self.makeDataRequest(_filtered_smiles) # make call for data!
			_measured_data = json.loads(_response.content)
			# _measured_data.update(json.loads(_response.content))
			logging.info("Measured Data: {}".format(_measured_data))
		except Exception as e:
			logging.warning("Exception making request to Measured: {}".format(e))
			_response_dict.update({'data': "data not found"})
			for prop in request_dict['props']:
				_response_dict.update({'prop': prop})
				self.redis_conn.publish(request_dict['sessionid'], json.dumps(_response_dict))
			return

		logging.info("Measured Data: {}".format(_measured_data))
		logging.info("Request props: {}".format(request_dict['props']))
			

		# get requested properties from results:
		for prop in request_dict['props']:
			try:
				_data_obj = self.getPropertyValue(prop, _measured_data)
				logging.info("data object: {}".format(_data_obj))
				_response_dict.update(_data_obj)
				
				# push one result at a time if node/redis:
				_result_json = json.dumps(_response_dict)
				self.redis_conn.publish(request_dict['sessionid'], _result_json)

			except Exception as err:
				logging.warning("Exception occurred getting Measured data: {}".format(err))
				_response_dict.update({'data': "cannot reach {} calculator".format(request_dict['calc'])})

				logging.info("##### session id: {}".format(request_dict['sessionid']))

				# node/redis stuff:
				self.redis_conn.publish(request_dict['sessionid'], json.dumps(_response_dict))