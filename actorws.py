import requests
import json
import logging



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
		self.dsstox_result_keys = ['casrn', 'dsstoxSubstanceId', 'preferredName', 'smiles', 'iupac']
		self.result_obj = {
			'calc': "actorws",
			'prop': "",
			'data': {},
		}


	def make_request(self, url, payload):
		try:
			_response = requests.get(url, params=payload, timeout=10)
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
			return {'success': False, 'error': "error connecting to actorws", 'data': None} 
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

		_dsstox_results = self.make_request(self.dsstox_url, _payload)

		logging.warning("DSSTOX RESULTS: {}".format(_dsstox_results))

		try:
			_dsstox_results = _dsstox_results['DataList']['list'][0]
		except Exception as e:
			logging.warning("Error getting dsstox results key:vals: {}".format(e))
			logging.warning("Using only Jchem WS results instead..")
			# raise e  # raise it for now

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
		_chemid_results = self.make_request(self.chemid_url, {'identifier': chemical})

		logging.warning("CHEMID RESULTS: {}".format(_chemid_results))

		try:
			_chemid_results = _chemid_results['DataRow']
		except KeyError as e:
			logging.warning("'DataRow' key not found in chemid results.. Returning None..")
			return None

		_results = self.result_obj
		_results['prop'] = "chemid"
		_result_key = self.chemid_result_keys[0]  # only one key needed for now
		
		if _result_key in _chemid_results:
			_results['data'].update({'gsid': _chemid_results.get(_result_key)})  # getting synGsid key:val

		# todo: add more error handling, waiting to start new cheminfo workflow w/ actorws first..

		return _results

	##########################################################################