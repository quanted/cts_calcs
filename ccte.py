import requests
import logging
import os
import html
import json


PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))


class CCTE:
	"""
	Handles requests and responses for new public CCTE endpoints.
	Docs: https://api-ccte.epa.gov/docs
	"""

	def __init__(self):

		# Base URL and endpoints for public CCTE:
		self.ccte_base_url = "https://api-ccte.epa.gov/"
		self.chem_search_equal_url = "chemical/search/equal/{}"
		self.chem_search_substring_url = "chemical/search/contain/{}"
		self.chem_search_starting_url = "chemical/search/start-with/{}"
		self.chem_prop_dtxsid_url = "chemical/property/search/by-dtxsid/{}"
		self.chem_details_dtxsid_url = "chemical/detail/search/by-dtxsid/{}"
		self.chem_details_dtxcid_url = "chemical/detail/search/by-dtxcid/{}"
		self.chem_fate_url = "chemical/fate/search/by-dtxsid/{}"

		# Required headers for their API:
		self.headers = {
			"Accept": "application/json",
			"x-api-key": None  # API key they provided us (set in env/config)
		}

		self.dsstox_result_keys = ['casrn', 'dtxsid', 'preferredName', 'smiles']  # result keys for chem info
		self.keys_of_interest = ["dtxsid", "casrn", "preferredName", "smiles"]  # (from actorws.py)

		# Data retrieved if chemaxon isn't available:
		self.chemid_keys_map = {
			'casrn': 'casrn',
			'preferredName': 'preferredName',
			'synGsid': 'gsid',
			'dtxsid': 'dsstoxSubstanceId',
			'dtxcid': 'dtxcid',
			'smiles': 'smiles',
			'molFormula': 'formula',
			'averageMass': 'mass',
			'monoisotopicMass': 'exactMass',
			'iupacName': 'iupac'
		}

		self.models_path = os.path.join(PROJECT_ROOT, "models")

		# Map for CCTE endpoints:
		self.chem_result_map = {
			"search": {
				"filename": "ccte_public_chemical_search_results_example.json",
				"data": None,
				"keys": ["dtxsid", "casrn", "preferredName", "smiles"]  # keys for cts
			},
			"details": {
				"filename": "ccte_public_chemical_details_results_example.json",
				"data": None
			},
			"fate": {
				"filename": "ccte_public_chemical_fate_results_example.json",
				"data": None
			},
			"property": {
				"filename": "ccte_public_chemical_property_results_example.json",
				"data": None
			}
		}

		# Wrapped result for CTS:
		self.wrapped_results = {
			"calc": None,
			"prop": None,
			"data": None
		}

		self.set_result_keys()

		self.set_api_key()

	def set_result_keys(self):
		"""
		Reads in CCTE json files of example results and adds
		them to the result map.
		"""
		# TODO: Add a "search" model for results as well.
		for key, result_obj in self.chem_result_map.items():
			filename = result_obj.get("filename")
			with open(os.path.join(self.models_path, filename)) as f:
				result_obj["data"] = json.load(f)

	def set_api_key(self):
		"""
		Sets API key from environment/config for request header.
		"""
		api_key = os.getenv("CCTE_API_KEY") 
		if not api_key:
			raise Exception("CCTE_API_KEY not found in environment.\n\
				Needed for making requests to public CCTE server (https://api-ccte.epa.gov/docs)")
		self.headers["x-api-key"] = api_key

	def wrap_results(self, results):
		"""
		Wraps validated and curated response data from CCTE for
		use by CTS (calc, prop, data keys, etc.).
		"""
		cts_results = dict(self.wrapped_results)
		cts_results["calc"] = "actorws"  # (from actorws.py)
		cts_results["prop"] = "dsstox"  # (from actorws.py)
		cts_results["data"] = results
		return cts_results

	def validate_response(self, response):
		"""
		Validates response from CCTE endpoints.
		Returns resonse content as object.
		"""
		data_obj = None
		if response.status_code != 200:
			logging.warning("ccte request non-200: {}".format(response))
			return {"status": False, "error": "Error making request to CCTE."}
		try:
			data_obj = json.loads(response.content)
		except Exception as e:
			logging.error("Cannot serialize response content: {}".format(e))
			return {"status": False, "error": "Error parsing CCTE data."}
		return data_obj

	def make_search_request(self, chemical):
		"""
		Makes chemical search request for getting chemical info.
		Accepts DTXSID, DTXCID, InChlKey, name, and CASRN as inputs.
		"""
		url = self.ccte_base_url + self.chem_search_equal_url.format(html.escape(chemical))
		response = None
		try:
			response = requests.get(url, headers=self.headers)
		except Exception as e:
			logging.warning("ccte make_search_request exception, url: {}: {}".format(url, e))
			return False
		response_obj = self.validate_response(response)
		if not isinstance(response_obj, list) and response_obj.get("status") != True:
			# TODO: More exception handling?
			return response_obj
		if len(response_obj) != 1:
			logging.warning("More than one chemical returned in chemical search: {}".format(response_obj))
			# TODO: idk, gonna pick the first one for now.
		results = self.get_search_results(response_obj)
		return self.wrap_results(results)

	def get_search_results(self, response_obj):
		"""
		Gets the keys CTS needs from CCTE chemical search response.
		"""
		data_obj = response_obj[0]  # expecting single-item in list for search results
		data_for_cts = {key: val for key, val in data_obj.items() if key in self.chem_result_map["search"]["keys"]}
		return data_for_cts

	def make_details_request(self, dtxsid):
		"""
		Makes chemical details request. Needs a DTXSID as
		the input.
		"""
		try:
			url = self.ccte_base_url + self.chem_details_dtxsid_url.format(html.escape(dtxsid))
			response = requests.get(url, headers=self.headers)
			logging.info("Search response for {}: {}".format())
			return json.loads(response.content)
		except Exception as e:
			logging.warning("ccte make_details_request exception, url: {}: {}".format(url, e))
			return False

	def get_details_results(self, response):
		"""
		Gets chemical details results.
		"""
		pass

	def make_propery_request(self, dtxsid):
		"""
		Makes a chemical property request using DTXSID
		as the input type.
		"""
		try:
			url = self.ccte_base_url + self.chem_prop_dtxsid_url.format(html.escape(dtxsid))
			response = requests.get(url, headers=self.headers)
			logging.info("Search response for {}: {}".format())
			return json.loads(response.content)
		except Exception as e:
			logging.warning("ccte make_propery_request exception, url: {}: {}".format(url, e))
			return False

	def get_property_results(self, response):
		"""
		Gets chemical property results.
		"""
		pass

	def make_fate_request(self, dtxsid):
		"""
		Makes a chemical search request using DTXSID as the
		input type.
		"""
		try:
			url = self.ccte_base_url + self.chem_fate_url.format(html.escape(dtxsid))
			response = requests.get(url, headers=self.headers)
			logging.info("Search response for {}: {}".format())
			return json.loads(response.content)
		except Exception as e:
			logging.warning("ccte make_fate_request exception, url: {}: {}".format(url, e))
			return False

	def get_fate_results(self, response):
		"""
		Gets chemical fate results.
		"""
		pass
