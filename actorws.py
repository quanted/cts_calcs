import requests
import json
import logging



class CCTE:
	"""
	New endpoints to use instead of ACTORWS.
	https://ccte-api-ccd.epa.gov
	"""
	def __init__(self):
		self.base_url = "https://ccte-api-ccd.epa.gov"
		self.batch_url = self.base_url + "/batchsearch/results"
		self.single_url = self.base_url + "/ccdapp1/search/chemical"
		self.calc = "actorws"  # NOTE: Test before changing this
		self.props = ['dsstox', 'chemid']
		self.batch_result_keys = ["input", "foundBy", "dtxsid", "dtxcid", "casrn", "preferredName"]
		self.single_result_keys = ["dtxsid", "dtxcid", "searchMatch", "rank", "hasStructureImage", "searchWord"]

		self.keys_of_interest = ["dtxsid", "casrn", "preferredName"]
		# self.dsstox_result_keys = ['casrn', 'dsstoxSubstanceId', 'preferredName', 'smiles', 'iupac']  # what we used to get back

		self.batch_request = {
			"inputType": "IDENTIFIER",
			"identifierTypes": [],
			"searchItems": ""
		}

		########################################
		# BATCH REQUEST EXAMPLES
		########################################
		# Example batch requests 1:
		# {
		#     "inputType" : "IDENTIFIER",
		# 	"identifierTypes" : ["CHEMICAL_NAME"],
		# 	"searchItems" : "TOLMETIN SODIUM\nSODIUM SUCCINATE\nSODIUM SUCCINATE\nSELENIUM SULFIDE\nbpa"
		# }
		# Example batch requests 2:
		# {
		#     "inputType" : "IDENTIFIER",
		# 	"identifierTypes" : ["CASRN"],
		# 	"searchItems" : "3/1/4860\n75-O5-8\n11121-31-6\n7782-50-5;\n0000542881\n0000542881;\n75058\n75-05-8\n75-50-8\n75-08-5"
		# }
		# Example batch requests 3:
		# {
		#     "inputType" : "IDENTIFIER",
		# 	"identifierTypes" : ["DTXSID"],
		# 	"searchItems" : "DTXSID8021569\nDTXSID8021640\nDTXSID8021642\nDTXSID8021644\nDTXSID8021646\nDTXSID8021690\nDTXSID8023923\nDTXSID8039241\nDTXSID8044466\nDTXSID8060955\nDTXSID8075049"
		# }
		# Example batch requests 4:
		# {
		# 	"inputType" : "MASS",
		#     "identifierTypes" : [""],
		# 	"searchItems" : "189\n191",
		#     "massError" : 0.5
		# }
		# Example batch requests 5:
		# {
		# 	"inputType" : "EXACT_FORMULA",
		#     "identifierTypes" : [""],
		# 	"searchItems" : "XYX\nC6H6O"
		# }
		# Example batch requests 6:
		# {
		# 	"inputType" : "MSREADY_FORMULA",
		#     "identifierTypes" : [""],
		# 	"searchItems" : "XYX\nC6H6O"
		# }
		########################################

		########################################
		# SINGLE REQUEST EXAMPLES
		########################################
		# Single request example 1:
		# 	https://ccte-api-ccd.epa.gov/ccdapp1/search/chemical/start-with/benzene
		# Single request example 2:
		# 	https://ccte-api-ccd.epa.gov/ccdapp1/search/chemical/equal/DTXSID0020103
		# Single request example 3:
		# 	https://ccte-api-ccd.epa.gov/ccdapp1/search/chemical/equal/120155-79-5
		# Single request example 4:
		# 	https://ccte-api-ccd.epa.gov/ccdapp1/search/chemical/contain/citrate
		########################################

		self.result_obj = {
			"calc": "actorws",
			"prop": "",
			"data": {},
		}

	def _make_request(self, url, data):
		try:
			_response = requests.post(url, data=data, headers={'Content-Type': 'application/json'}, timeout=15)
		except requests.exceptions.Timeout as e:
			logging.warning("Request to {} timed out.. No data from actorws..".format(url))
			return None
		except requests.exceptions.ConnectionError as e:
			logging.warning("Connection error for {}.. No data from actorws..".format(url))
			return None
		except Exception as e:
			logging.warning("Exception occurred in chemical information module: {}".format(e))
			return None
			
		if _response.status_code != 200:
			# return {'success': False, 'error': "error connecting to actorws", 'data': None}
			logging.warning("Exception in actorws.py making request to actorws: {}".format(_response))
			raise Exception("ACTORWS request was not successful.")

		return json.loads(_response.content)

	def get_chemical_results(self, chemical, chem_type):
		"""
		Makes request to CCTE batch endpoint to get chemical data.
		
		NOTE: Worked with chemical name or CASRN as input, not so much SMILES.
		TODO: Use chemical name or CAS from chemaxon if user-entered chemical is not
		one of those two.
		"""

		# Builds request object:
		post_data = dict(self.batch_request)
		if chem_type == "name":
			post_data["identifierTypes"] = ["CHEMICAL_NAME"]
		elif chem_type == "casrn":
			post_data["identifierTypes"] = ["CASRN"]
		else:
			return  # TODO: Handle excpetion
		post_data["searchItems"] = chemical

		# Makes request to CCTE:
		try:
			results = self._make_request(self.batch_url, post_data)
			results = results[0]
		except Exception as e:
			logging.warning("Error getting CCTE data: {}".format(e))
			logging.warning("Using only Jchem WS results instead, setting dsstox value to N/A")
			_results = dict(self.result_obj)
			_results['prop'] = "dsstox"
			_results['data'] = {'dsstox': "N/A"}
			return _results

		# Gets results of interest:
		_results = dict(self.result_obj)
		_results['prop'] = "dsstox"
		for key, val in results.items():
			if key in self.keys_of_interest:
				_results['data'][key] = val

		return _results



class ACTORWS(object):
	"""
	Uses actor web services to obtain curated
	CAS#, SMILES, preferred name, iupac, and DTXSID.
	Location: https://actorws.epa.gov/actorws/
	"""
	def __init__(self):
		self.base_url = "https://actorws.epa.gov/actorws"
		self.chemid_url = "https://actorws.epa.gov/actorws/chemIdentifier/v01/resolve.json"  # ?identifier=[chemical name or SMILES]
		self.dsstox_url = "https://actorws.epa.gov/actorws/dsstox/v02/casTable.json"  # ?identifier=[casrn or gsid]
		self.calc = "actorws"
		self.props = ['dsstox', 'chemid']
		self.chemid_result_keys = ['synGsid']  # chemidentifier result key of interest
		self.chemid_all_keys = ['origIdentifier', 'casrn', 'preferredName', 'synGsid', 'synType', 'synIdentifier', 'dtxsid', 'dtxcid', 'jChemInChIKey', 'indigoInChIKey', 'smiles', 'molFormula', 'molWeight', 'collidingGsid', 'collidingCasrn', 'collidingPreferredName', 'trimmedWhitespace', 'trimmedLeadingZeros', 'reformattedIdentifier', 'checksum', 'processedAs', 'infoMsg', 'warningMsg', 'msReadyForms', 'qsarForms', 'imageURL']
		
		self.chemid_keys_map = {
			'casrn': 'casrn',
			'preferredName': 'preferredName',
			'synGsid': 'gsid',
			'dtxsid': 'dsstoxSubstanceId',
			'dtxcid': 'dtxcid',
			'smiles': 'smiles',
			'molFormula': 'formula',
			'molWeight': 'mass',
		}

		self.dsstox_result_keys = ['casrn', 'dsstoxSubstanceId', 'preferredName', 'smiles', 'iupac']
		self.result_obj = {
			'calc': "actorws",
			'prop': "",
			'data': {},
		}


	def make_request(self, url, payload):
		try:
			_response = requests.get(url, params=payload, timeout=15)
		except requests.exceptions.Timeout as e:
			logging.warning("Request to {} timed out.. No data from actorws..".format(url))
			return None
		except requests.exceptions.ConnectionError as e:
			logging.warning("Connection error for {}.. No data from actorws..".format(url))
			return None
		except Exception as e:
			logging.warning("Exception occurred in chemical information module: {}".format(e))
			return None
			
		if _response.status_code != 200:
			# return {'success': False, 'error': "error connecting to actorws", 'data': None}
			logging.warning("Exception in actorws.py making request to actorws: {}".format(_response))
			raise Exception("ACTORWS request was not successful.")

		return json.loads(_response.content)



	############### "PUBLIC" METHODS BELOW ########################

	def get_dsstox_results(self, chemical, id_type):
		"""
		Makes request to actowws dsstox for the following
		result keys: casrn, dsstoxSubstanceId, preferredName, smiles, and iupac
		Input: cas number or gsid obtained from actorws chemicalIdentifier endpoint
		Output: Dictionary of above result key:vals
		"""
		_payload = {}
		if id_type == 'gsid':
			_payload = {'gsid': chemical}
		elif id_type == 'CAS#':
			_payload = {'casrn': chemical}

		try:
			_dsstox_results = self.make_request(self.dsstox_url, _payload)
			_dsstox_results = _dsstox_results['DataList']['list'][0]
		except Exception as e:
			logging.warning("Error getting dsstox results key:vals: {}".format(e))
			logging.warning("Using only Jchem WS results instead, setting dsstox value to N/A")
			_results = self.result_obj
			_results['prop'] = "dsstox"
			_results['data'] = {'dsstox': "N/A"}
			return _results

		_results = self.result_obj
		_results['prop'] = "dsstox"
		for _key, _val in _dsstox_results.items():
			if _key in self.dsstox_result_keys:
				_results['data'][_key] = _val

		return _results


	def get_chemid_results(self, chemical):
		"""
		Makes request to actorws chemicalIdentifier endpoint for
		'synGsid' to be used for dsstox if cas isn't provided by user.

		Inputs: chemical - either a chemical name or smiles
		Output: Dictionary with results_obj keys and synGsid
		"""
		try:
			_chemid_results = self.make_request(self.chemid_url, {'identifier': chemical})
			_chemid_results = _chemid_results['DataRow']
		except Exception as e:
			logging.warning("Exception getting chemid results from actorws: {}".format(e))
			# return None
			return {}

		return _chemid_results

	##########################################################################