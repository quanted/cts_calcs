import requests
import json
import logging
import os
import re
from collections import defaultdict

from .calculator import Calculator
from .chemical_information import SMILESFilter
from .ccte import CCTE
from .mongodb_handler import MongoDBHandler



headers = {'Content-Type': 'application/json'}

db_handler = MongoDBHandler()



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

		self.ccte_fate_map = {
			"Soil Adsorp. Coeff. (Koc)": "koc",
			"Bioconcentration Factor": "log_bcf",
			"Bioaccumulation Factor": "log_baf"  # NOTE: unverified key name
		}

		self.fate_result = {
			"valueType": "experimental",
			"endpointName": None
		}

		self.props = ["melting_point", "boiling_point", "water_sol", "vapor_press", "henrys_law_con", "kow_no_ph"]
		self.fate = ["koc", "log_bcf", "log_baf"]

		self.result_structure = {
			'structure': '',
			'propertyname': '',
			'propertyvalue': None
		}

		self.response_obj = {
			'calc': "measured",  # todo: change to metabolizer, change in template too
			'prop': "pchem",
			'data': None,
			'chemical': None,
			'request_post': None            
		}

	def getPostData(self):
		return {"structure": ""}


	def add_cts_keys_prop_data(self, results):
		"""
		Curates prop data from CCTE into keys that CTS understands,
		e.g., 'prop', 'data', 'method', 'chemical'.
		"""
		new_results = []
		for data_obj in results:
			if not data_obj["propertyId"] in list(self.ccte_prop_map.keys()):
				continue
			new_data_obj = dict(data_obj)
			new_data_obj["prop"] = self.ccte_prop_map[data_obj["propertyId"]]
			new_data_obj["method"] = self.convert_to_acronym(data_obj["source"])
			new_data_obj["data"] = data_obj["value"]
			new_results.append(new_data_obj)
		return new_results


	def add_cts_keys_fate_data(self, results):
		"""
		Curates prop data from CCTE into keys that CTS understands,
		e.g., 'prop', 'data', 'method', 'chemical'.
		"""
		new_results = []
		for data_obj in results:
			if not data_obj["valueType"] == "experimental" \
				or not data_obj["endpointName"] in list(self.ccte_fate_map.keys()):
					continue
			new_data_obj = dict(data_obj)
			new_data_obj["prop"] = self.ccte_fate_map[data_obj["endpointName"]]
			new_data_obj["method"] = self.convert_to_acronym(data_obj["modelSource"])
			new_data_obj["data"] = data_obj["resultValue"]
			new_results.append(new_data_obj)
		return new_results


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


	def parse_pka_data(self, db_results):
		"""
		Parses db results into response for speciation workflow.
		"""
		data_obj = {}
		pka_list = []

		for i in range(1, 7):
			key = "pKa_{}".format(i)
			val = db_results[key]
			if val == "":
				break
			pka_list.append(val)
		
		smiles = db_results.get("Standardized_SMILES")
		if not smiles:
			smiles = db_results.get("SMILES")
		
		data_obj["pka_list"] = pka_list
		data_obj["smiles"] = smiles
		data_obj["status"] = True
		
		return data_obj


	def handle_pka_request(self, request_dict):
		"""
		Requests Measured pka data from mongodb. Uses "DTXSID" to find data
		or "Standardized_SMILES" column if the former doesn't exist.

		Example DB document:
		{
			'_id': ObjectId('671fa9b3fdf54bb02fa35eaa'),
			'SheetName': 'Prankerd',
			'Name': 'Aspirin',
			'CASRN': '50-78-2',
			'CASRN_updated': '50-78-2',
			'DTXSID': 'DTXSID5020108',
			'Standardized_SMILES': '',
			'SMILES': 'CC(=O)OC1=C(C=CC=C1)C(O)=O',
			'Temp °C': 17,
			'pKa_1': 3.565,
			'pKa_2': '',
			'pKa_3': '',
			'pKa_4': '',
			'pKa_5': '',
			'pKa_6': '',
			'Compound_Type': '',
			'Condition': '',
			'Ionization': '-',
			'Reference': 'Edwards LJ, The hydrolysis of aspirin, Trans. Farad. Soc., 46, 723–735\n(1950).',
			'Ref_number': ''
		}
		"""

		dtxsid = request_dict.get("dtxsid")

		smiles = request_dict.get("chemical")  # TODO: Determine best key to use, may be "smiles"

		try:

			db_handler.connect_to_db()

			# try:
			
			if not db_handler.is_connected:
				logging.warning("OPERA DB not connected.")
				return False

			# TODO: Check that dtxsid value already in request_dict.
			# NOTE: Use standardized smiles if dtxsid not available in DB.

			db_results = None

			# db_results = db_handler.find_pka_document
			db_results = db_handler.pka_collection.find_one({
				"DTXSID": dtxsid
			})

			if len(db_results) < 1:
				# Checks to see if Standardized_SMILES exists if no results from dtxsid:
				db_results = db_handler.pka_collection.find_one({
					"Standardized_SMILES": smiles
				})

			if len(db_results) > 0:
				db_results = self.parse_pka_data(db_results)
				if db_results and "_id" in db_results:
					del db_results["_id"]
				return db_results
			else:
				return False

		except Exception as e:
			logging.error("calculator_measured handle_pka_request error: {}".format(e))
			return False

		finally:
			db_handler.mongodb_conn.close()

		return db_results


	def data_request_handler(self, request_dict):

		# logging.warning("calculator_measured request_dict: {}".format(request_dict))

		_filtered_smiles = ''
		_response_dict = {}
		_measured_data = {}

		# fill any overlapping keys from request:
		for key in request_dict.keys():
			if not key == 'nodes':
				_response_dict[key] = request_dict.get(key)
		_response_dict.update({'request_post': request_dict, 'method': None})

		# NOTE: Assuming no SMILES filter by calc for pka in speciation:
		if request_dict.get("prop") == "ion_con" or request_dict.get("service") == "getSpeciationData":
			_response_obj = dict(self.response_obj)
			try:
				# Gets pka from measured DB
				db_results = self.handle_pka_request(request_dict)
				_response_obj['data'] = db_results
				_response_obj['chemical'] = request_dict.get('chemical')
				_response_obj['request_post'] = request_dict
			except Exception as e:
				logging.error("calculator_measured exception: {}".format(e))
				_response_obj.update({"valid": False, 'error': "Error getting data from measured"})
				return _response_obj

			return _response_obj

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



		####################################################################
		# TODO: Add conditional to only request props or fate endpoints if
		# they're in the user request.
		####################################################################

		prop_response = None
		_response_dict.update({"prop_results": []})

		if any(prop in request_props for prop in self.props):
			# Makes property request to CCTE for MP, BP, WS, VP, HL, and KOW.
			prop_response = self.make_propery_request(dtxsid)
			if not prop_response:
				_response_dict.update({
					'data': "Cannot get prop from CCTE",
					'valid': False
				})
				return _response_dict
			prop_results = self.get_property_results(prop_response)
			curated_results = self.add_cts_keys_prop_data(prop_results)
			final_prop_results = self.group_by_acronym(curated_results)
			_response_dict["prop_results"] = _response_dict["prop_results"] + final_prop_results

		if any(fate in request_props for fate in self.fate):
			# Makes fate request to CCTE for KOC, BCF, and BAF.
			fate_response = self.make_fate_request(dtxsid)
			if not fate_response:
				_response_dict.update({
					'data': "Cannot get fate data from CCTE",
					'valid': False
				})
				return _response_dict

			fate_results = self.add_cts_keys_fate_data(fate_response)
			_response_dict["prop_results"] = _response_dict["prop_results"] + fate_results

		return _response_dict
