import requests
import json
import logging
import os
import re
from collections import defaultdict

from .calculator import Calculator
from .chemical_information import SMILESFilter
from .ccte import CCTE


headers = {'Content-Type': 'application/json'}



class MeasuredCalc(Calculator, CCTE):
	"""
	Measured Calculator 2
	Gets values from CCTE endpoints.
	"""

	def __init__(self):
		Calculator.__init__(self)

		CCTE.__init__(self)

		self.postData = {"smiles" : ""}
		self.name = "measured"
		self.request_timeout = 20
		self.melting_point = 0.0

		# map workflow parameters to test
		self.propMap = {
			'melting_point': {
			   'result_key': 'melting_point',
			   'prop_id': 'melting-point'
			},
			'boiling_point': {
			   'result_key': 'boiling_point',
			   'prop_id': 'boiling-point'
			},
			'water_sol': {
			   'result_key': 'water_solubility',
			   'prop_id': 'water-solubility'
			},
			'vapor_press': {
			   'result_key': 'vapor_pressure',
			   'prop_id': 'vapor-pressure'
			},
			'henrys_law_con': {
				'result_key': 'henrys_law_constant',
				'prop_id': 'henrys-law'
			},
			'kow_no_ph': {
				'result_key': 'log_kow',
				'prop_id': ''
			},
			'koc': {
				'result_key': 'log_koc',
				'prop_id': 'log-koc',
				'endpointName': "Soil Adsorp. Coeff. (Koc)"
			},
			'log_bcf': {
				'result_key': 'log_bcf',
				'prop_id': 'log-bcf',
				'endpointName': "Bioconcentration Factor"
			},
			'log_baf': {
				'result_key': 'log_baf',
				'prop_id': 'log-baf',
				'endpointName': "Bioaccumulation Factor"	
			}
		}

		# ccte propertyId to CTS prop map:
		self.ccte_prop_map = {
			"melting-point": "melting_point",
			"boiling-point": "boiling_point",
			"water-solubility": "water_sol",
			"vapor-pressure": "vapor_press",
			"henrys-law": "henrys_law_con",
			"logkow-octanol-water": "kow_no_ph"
		}

		self.props = ["melting_point", "boiling_point", "water_sol", "vapor_press", "henrys_law_con", "kow_no_ph"]
		self.fate = ["koc", "log_bcf", "log_baf"]

		self.result_structure = {
			'structure': '',
			'propertyname': '',
			'propertyvalue': None
		}

	def getPostData(self):
		return {"structure": ""}


	def add_cts_keys(self, results):
		"""
		Curates prop data from CCTE into keys that CTS understands,
		e.g., 'prop', 'data', 'method', 'chemical'.
		"""
		for data_obj in results:
			if not data_obj["propertyId"] in list(self.ccte_prop_map.keys()):
				continue
			data_obj["prop"] = self.ccte_prop_map[data_obj["propertyId"]]
			data_obj["method"] = self.convert_to_acronym(data_obj["source"])

			# TODO: Value conversions where necessary.

			data_obj["data"] = data_obj["value"]

		return results


	def group_by_acronym(self, props_list):
		"""
		Consolidates data objects with the same prop and method
		into one object and concatenates data into a comma-separated string.
		"""
		new_list = []
		new_dict = {}
		# Process each entry in the data list
		for item in props_list:

			if not item["propertyId"] in list(self.ccte_prop_map.keys()):
				continue

			match_key = (item['method'], item['prop'])
			if match_key not in new_dict:
				# Initializes dict object if it doesn't exist already:
				new_dict[match_key] = {k: v for k, v in item.items() if k != "data"}
				new_dict[match_key]["data"] = str(item["data"])
			else:
				# Concatenates "data" in comma-separated string:
				new_dict[match_key]["data"] += f", {item['data']}"
		new_list = list(new_dict.values())
		return new_list


	def validate_response(self, response):
		"""
		Validates sparc response.
		Returns False if data is null, or any other
		values that indicate an error
		"""
		if response.status_code != 200:
			logging.warning("measured server response status: {}".format(response.status_code))
			return False

		# successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
		# json_obj = json.loads(response.content)

		# TODO: verify if blank data, finding the source of the empty water sol values...
		return True


	def convert_to_acronym(self, method):

		ignore_words = ["et", "et.", "al", "al.", "and", "of", "the", "at"]

		# Remove text inside parentheses
		cleaned_method = re.sub(r'[^A-Za-z]', '', method).strip()
		
		# Check if it's a camel case or space-separated string
		if ' ' in cleaned_method:
			words = [word for word in cleaned_method.split() if word.lower() not in ignore_words]
			acronym = ''.join(word[0].upper() for word in words if word[0].isalpha())
		else:
			# CamelCase or mixed case: Extract only the uppercase letters
			acronym = ''.join(char for char in cleaned_method if char.isupper() and char.isalpha())
		return acronym


	def data_request_handler(self, request_dict):

		_filtered_smiles = ''
		_response_dict = {}
		_measured_data = {}

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

		chem_info = request_dict.get("chem_info", {})
		dtxsid = chem_info.get("dtxsid")
		pchem_request = request_dict.get("pchem_request")
		request_props = pchem_request.get("measured", {})

		if not request_props:
			logging.warning("calculator_measured no request props: {}".format(pchem_request))
			_response_dict.update({
				'data': "Cannot get properties from CCTE",
				'valid': False
			})
			return _response_dict

		# Makes property request to CCTE for MP, BP, WS, VP, HL, and KOW.
		prop_response = self.make_propery_request(dtxsid)

		if not prop_response:
			_response_dict.update({
				'data': "Cannot get properties from CCTE",
				'valid': False
			})
			return _response_dict

		prop_results = self.get_property_results(prop_response)
		full_results = self.add_cts_keys(prop_results)
		results = self.group_by_acronym(full_results)
		_response_dict.update({"prop_results": results})


		# TODO: 
		# # 2. Make fate request to CCTE for KOC, BCF, and BAF.
		# fate_response = self.make_fate_request(dtxsid)
		# logging.warning("calculator_measured fate_response: {}".format(fate_response))
		# if not fate_response:
		# 	_response_dict.update({
		# 		'data': "Cannot get fate data from CCTE",
		# 		'valid': False
		# 	})
		# 	return _response_dict


		return _response_dict
