__author__ = "np"

import requests
import logging
import json
from .calculator import Calculator
from .jchem_properties import Tautomerization


class ChemInfo(object):
	"""
	Suggested class for organizing chemical info in CTS.
	Captures the objects and key:vals used for obtaining
	chemical data.
	"""
	def __init__(self, chemical=None):
		self.chem_obj = {
			'chemical': chemical,  # user-entered molecule of any chemaxon format
			'orig_smiles': "",  # 'chemical' converted to smiles (pre-filtering)
			'smiles': "",  # result of orig_smiles after cts filtering (standardized smiles)
			'preferredName': "",
			'iupac': "",
			'formula': "",
			'cas': "",
			'dtxsid': "",
			'mass': "",
			'exactMass': ""
		}

	def getChemInfo(self, request_post):
		"""
		Initial attempt for a single chem info function that's used
		by tasks.py in cts celery and cts_rest.py in cts api

		NOTE/TODO: Currently doesn't use the actorws services, like
		cts_rest.py does for the Chemical Editor. Use func from cts_rest
		once this works..
		"""
		chemical = request_post.get('chemical')
		get_sd = request_post.get('get_structure_data')  # bool for getting <cml> format image for marvin sketch

		response = self.convertToSMILES({'chemical': chemical})
		orig_smiles = response['structure']
		filtered_smiles_response = smilesfilter.filterSMILES(orig_smiles)
		filtered_smiles = filtered_smiles_response['results'][-1]

		logging.info("Filtered SMILES: {}".format(filtered_smiles))

		jchem_response = self.getChemDetails({'chemical': filtered_smiles})  # get chemical details

		molecule_obj = {'chemical': filtered_smiles}
		for key, val in jchem_response['data'][0].items():
			molecule_obj[key] = val
			# chem_list.append(molecule_obj)

		if request_post.get('is_node'):
			#### only get these if gentrans single mode: ####
			molecule_obj.update({'node_image': self.nodeWrapper(filtered_smiles, MetabolizerCalc().tree_image_height, MetabolizerCalc().tree_image_width, MetabolizerCalc().image_scale, MetabolizerCalc().metID,'svg', True)})
			molecule_obj.update({
				'popup_image': self.popupBuilder(
					{"smiles": filtered_smiles}, 
					MetabolizerCalc().metabolite_keys, 
					"{}".format(request_post.get('id')),
					"Metabolite Information")
			})
			##################################################

		wrapped_post = {
			'status': True,  # 'metadata': '',
			'data': molecule_obj,
			'request_post': request_post
		}

		logging.info("Returning Chemical Info: {}".format(json_data))

		return wrapped_post  # return json object!




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
        _response = requests.get(url, params=payload, timeout=10)
        if _response.status_code != 200:
            return {'success': False, 'error': "error connecting to actorws", 'data': None} 
        return json.loads(_response.content)


    ##### "PUBLIC" METHODS BELOW #####
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
        except KeyError as e:
            logging.warning("Error getting dsstox results key:vals..")
            raise e  # raise it for now

        # what key:vals should be with results??

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
            logging.warning("Error getting dsstox results key:vals..")
            raise e  # raise it for now

        # what key:vals should be with results??

        _results = self.result_obj
        _results['prop'] = "chemid"
        _result_key = self.chemid_result_keys[0]  # only one key needed for now
        
        if _result_key in _chemid_results:
            _results['data'].update({'gsid': _chemid_results.get(_result_key)})  # getting synGsid key:val

        # todo: add more error handling, waiting to start new cheminfo workflow w/ actorws first..

        return _results




class SMILESFilter(object):
	"""
	This is the smilesfilter.py module as a class and
	clumped in with other classes related to chem info.
	"""

	def __init__(self):
		self.max_weight = 1500  # max weight [g/mol] for epi, test, and sparc
		self.excludestring = {".","[Ag]","[Al]","[Au]","[As]","[As+","[B]","[B-]","[Br-]","[Ca]",
						"[Ca+","[Cl-]","[Co]","[Co+","[Fe]","[Fe+","[Hg]","[K]","[K+","[Li]",
						"[Li+","[Mg]","[Mg+","[Na]","[Na+","[Pb]","[Pb2+]","[Pb+","[Pt]",
						"[Sc]","[Si]","[Si+","[SiH]","[Sn]","[W]"}
		self.return_val = {
			"valid" : False,
			"smiles": "",
			"processedsmiles" : ""
		}

	def is_valid_smiles(self, smiles):

		_return_val = self.return_val
		_return_val['smiles'] = 

		if any(x in smiles for x in self.excludestring):
			return _return_val  # metal in smiles, invalid!

		try:
			processed_smiles = self.filterSMILES(smiles)
		except Exception as e:
			logging.warning("!!! Error in smilesfilter {} !!!".format(e))
			raise "smiles filter exception, possibly invalid smiles..."
				
		_return_val["valid"] = True
		_return_val["processedsmiles"] = processed_smiles

		return _return_val


	def singleFilter(self, request_obj):
		"""
		Calls single EFS Standardizer filter
		for filtering SMILES
		"""
		try:
			smiles = request_obj.get('smiles')
			action = request_obj.get('action')
		except Exception as e:
			logging.info("Exception retrieving mass from jchem: {}".format(e))
			raise
		post_data = {
			"structure": smiles,
			"actions": [
				action
			]
		}
		url = Calculator().efs_server_url + Calculator().efs_standardizer_endpoint
		return Calculator().web_call(url, post_data)


	def filterSMILES(self, smiles):
		"""
		cts ws call to jchem to perform various
		smiles processing before being sent to
		p-chem calculators
		"""

		# Updated approach (todo: more efficient to have CTSWS use major taut instead of canonical)
		# 1. CTSWS actions "removeExplicitH" and "transform".
		url = Calculator().efs_server_url + Calculator().efs_standardizer_endpoint
		post_data = {
			'structure': smiles,
			'actions': [
				"removeExplicitH",
				"transform"
			]
		}
		response = Calculator().web_call(url, post_data)

		filtered_smiles = response['results'][-1] # picks last item, format: [filter1 smiles, filter1 + filter2 smiles]
		
		# 2. Get major tautomer from jchem:
		taut_obj = Tautomerization()
		taut_obj.postData.update({'calculationType': 'MAJOR'})
		taut_obj.make_data_request(filtered_smiles, taut_obj)

		# todo: verify this is major taut result smiles, not original smiles for major taut request...
		major_taut_smiles = taut_obj.results['result']['structureData']['structure']

		# 3. Using major taut smiles for final "neutralize" filter:
		post_data = {
			'structure': major_taut_smiles, 
			'actions': [
				"neutralize"
			]
		}
		response = Calculator().web_call(url, post_data)

		final_smiles = response['results'][-1]
		logging.warning("FINAL FITERED SMILES: {}".format(final_smiles))

		return response


	def checkMass(self, chemical):
		"""
		returns true if chemical mass is less
		than 1500 g/mol
		"""
		logging.info("checking mass..")
		try:
			json_obj = Calculator().getMass({'chemical': chemical}) # get mass from jchem ws
		except Exception as e:
			logging.warning("!!! Error in checkMass() {} !!!".format(e))
			raise e
		logging.info("mass response data: {}".format(json_obj))
		struct_mass = json_obj['data'][0]['mass']
		logging.info("structure's mass: {}".format(struct_mass))

		if struct_mass < 1500  and struct_mass > 0:
			return True
		else:
			return False


	def clearStereos(self, smiles):
		"""
		clears stereoisomers from smiles
		"""
		try:
			response = self.singleFilter({'smiles':smiles, 'action': "clearStereo"})
			filtered_smiles = response['results'] # get stereoless structure
		except Exception as e:
			logging.warning("!!! Error in clearStereos() {} !!!".format(e))
			raise e
		return filtered_smiles


	def transformSMILES(self, smiles):
		"""
		N(=O)=O >> [N+](=O)[O-]
		"""
		try:
			response = self.singleFilter({'smiles':smiles, 'action': "transform"})
			filtered_smiles = response['results'] # get stereoless structure
		except Exception as e:
			logging.warning("!!! Error in transformSMILES() {} !!!".format(e))
			raise e
		return filtered_smiles


	def untransformSMILES(self, smiles):
		"""
		[N+](=O)[O-] >> N(=O)=O
		"""
		try:
			response = self.singleFilter({'smiles':smiles, 'action': "untransform"})
			filtered_smiles = response['results'] # get stereoless structure
		except Exception as e:
			logging.warning("!!! Error in untransformSMILES() {} !!!".format(e))
			raise e
		return filtered_smiles


	def parseSmilesByCalculator(self, structure, calculator):
		"""
		Calculator-dependent SMILES filtering!
		"""
		logging.info("Parsing SMILES by calculator..")
		filtered_smiles = structure

		#1. check structure mass..
		if calculator != 'chemaxon':
			logging.info("checking mass for: {}...".format(structure))
			if not self.checkMass(structure):
				logging.info("Structure too large, must be < 1500 g/mol..")
				raise "Structure too large, must be < 1500 g/mol.."

		#2-3. clear stereos from structure, untransform [N+](=O)[O-] >> N(=O)=O..
		if calculator == 'epi' or calculator == 'sparc' or calculator == 'measured':
			try:
				# clear stereoisomers:
				filtered_smiles = self.clearStereos(structure)
				logging.info("stereos cleared: {}".format(filtered_smiles))

				# transform structure:
				filtered_smiles = str(filtered_smiles[-1])
				filtered_smiles = str(self.untransformSMILES(filtered_smiles)[-1])
				logging.info("structure transformed..")
			except Exception as e:
				logging.warning("!!! Error in parseSmilesByCalculator() {} !!!".format(e))
				raise e

		# 4. Check for metals and stuff (square brackets):
		if calculator == 'epi' or calculator == 'measured':
			if '[' in filtered_smiles or ']' in filtered_smiles:
				# bubble up to calc for handling error
				raise Exception("{} cannot process metals...".format(calculator))

		return filtered_smiles