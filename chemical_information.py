__author__ = "np"

import requests
import logging
import json
from .calculator_metabolizer import MetabolizerCalc
from .actorws import ACTORWS, CCTE_EPA
from .smilesfilter import SMILESFilter
from .ccte import CCTE



class Molecule(object):
	"""
	Basic molecule object for CTS
	"""
	def __init__(self):

		# cts keys:
		self.chemical = ''  # initial structure from user (any chemaxon format)
		self.orig_smiles = ''  # before filtering, converted to smiles

		# chemaxon/jchem keys:
		self.smiles = ''  # post filtered smiles 
		self.formula = ''
		self.iupac = ''
		# self.cas = ''
		self.mass = ''
		self.structureData = ''
		self.exactMass = ''
		self.preferredName = ''	

	def createMolecule(self, chemical, orig_smiles, chem_details_response, get_structure_data=None):
		"""
		Gets Molecule attributes from Calculator's getChemDetails response
		"""
		# set attrs from jchem data:
		for key in self.__dict__.keys():
			if key != 'orig_smiles' and key != 'chemical':
				try:
					# if key == 'structureData' and get_structure_data == None:
					# 	pass
					if key == 'cas':
						# check if object with 'error' key instead of string of CAS#..
						if isinstance(chem_details_response['data'][0][key], dict):
							self.__setattr__(key, "N/A")
						else:
							self.__setattr__(key, chem_details_response['data'][0][key])	
					else:
						self.__setattr__(key, chem_details_response['data'][0][key])
				except KeyError as err:
					raise err
		# set cts attrs:
		self.__setattr__('chemical', chemical)
		self.__setattr__('orig_smiles', orig_smiles)

		return self.__dict__




class ChemInfo(object):
	"""
	Suggested class for organizing chemical info in CTS.
	Captures the objects and key:vals used for obtaining
	chemical data.
	"""
	def __init__(self, chemical=""):
		self.actorws_obj = ACTORWS()
		self.ccte_epa_obj = CCTE_EPA()
		self.ccte_obj = CCTE()
		self.smiles_filter_obj = SMILESFilter()
		self.calc_obj = MetabolizerCalc()  # note: inherits Calculator class as well
		self.cas_url = "https://cactus.nci.nih.gov/chemical/structure/{}/cas"  # associated CAS
		self.carbon_anomolies = {
			"C": "methane",
			"CC": "ethane",
			"CCC": "propane"
		}  # carbon-chain chems that comptox reads as a name instead of smiles (methane, ethane, and propane, respectively)
		self.chem_names_smiles_map = {
			"pfas": "Perfluorooctanoic acid"
		}
		self.chem_obj = [
			{
				'name': "chemical",
				'label': "Entered Chemical",
				'value': chemical
			},
			{
				'name': "orig_smiles",
				'label': "Initial SMILES",
				'value': ""
			},
			{
				'name': "smiles",
				'label': "Standardized SMILES",
				'value': ""
			},
			{
				'name': "preferredName",
				'label': "Preferred Name",
				'value': ""
			},
			{
				'name': "iupac",
				'label': "IUPAC",
				'value': ""
			},
			{
				'name': "formula",
				'label': "Formula",
				'value': ""
			},
			{
				'name': "casrn",
				'label': "Preferred CAS",
				'value': ""
			},
			{
				'name': "cas",
				'label': "Associated CAS",
				'value': ""
			},
			{
				'name': "dtxsid",
				'label': "DTXSID",
				'value': ""
			},
			{
				'name': "mass",
				'label': "Average Mass (g/mol)",
				'value': float
			},
			{
				'name': "exactmass",
				'label': "Monoisotopic Mass (g/mol)",
				'value': float
			}
		]
		self.wrapped_post = {
			'status': False,
			'data': None,
			'request_post': None
		}

	def create_cheminfo_table(self, workflow_obj):
		"""
		Creates object for workflow output page's user-input table.
		"""
		data = [
			{'Entered chemical': workflow_obj.chem_struct},
			{'Initial SMILES': workflow_obj.orig_smiles},
			{'Standardized SMILES': workflow_obj.smiles},
			{'IUPAC': workflow_obj.name},
			{'Formula': workflow_obj.formula},
			{'CAS #': workflow_obj.cas},
			{'Average Mass': workflow_obj.mass},
			{'Monoisotopic Mass': workflow_obj.exactMass}
		]
		return data

	def get_cheminfo(self, request_post, only_dsstox=False):
		"""
		Makes call to Calculator for chemaxon
		data. Converts incoming structure to smiles,
		then filters smiles, and then retrieves data
		:param request:
		:return: chemical details response json

		Note: Due to marvin sketch image data (<cml> image) being
		so large, a bool, "structureData", is used to determine
		whether or not to grab it. It's only needed in chem edit tab.
		"""
		# Updated cheminfo workflow with actorws:
		###########################################################################################
		# 1. Determine if user's chemical is smiles, cas, or drawn
		# 		a. If smiles, get gsid from actorws chemicalIdentifier endpoint
		#		b. If cas, get chem data from actorws dsstox endpoint
		#		c. If drawn, get smiles from chemaxon, then gsid like in a.
		# 2. Check if request from 1. "matched" (exist?)
		#		a. If 1b returns cas result, get cheminfo from dsstox results
		#		b. If 1a or 1c, use gsid from chemicalIdentifier and perform step 1b for dsstox cheminfo
		# 3. Use dsstox results: curated CAS#, SMILES, preferredName, iupac, and dtxsid
		#		a. Display in chemical editor.
		#############################################################################################

		chemical = request_post.get('chemical')
		get_sd = request_post.get('get_structure_data')  # bool for getting <cml> format image for marvin sketch
		is_node = request_post.get('is_node')  # bool for tree node or not
		_actor_results = {}  # dict for actorws results
		_gsid = None
		orig_smiles = None  # initial SMILES pre CTS filter
		is_name = False  # bool for whether smiles was actually acronym

		# Checks chemical against chem_name_smiles_map:
		chemical = self.check_name_smiles_map(chemical)

		is_valid_structure = self.check_structure_request(chemical)

		# Checks for valid structure:
		if "error" in is_valid_structure:
			response_obj = {}
			response_obj['status'] = False
			response_obj['error'] = is_valid_structure["error"]
			response_obj['request_post'] = request_post
			return response_obj

		# Determines chemical type from user (e.g., smiles, cas, name, etc.):
		chem_type = self.calc_obj.get_chemical_type(chemical)

		if not chem_type.get('type') or 'error' in chem_type:
			

			# TODO: Return error message if chem_type cannot be received from chemaxon.


			# If ChemAxon can't get chemical type, then try to get chem info
			# data from ACTORWS instead.
			logging.info("Couldn't get chemical type from ChemAxon. Trying to get data from ACTORWS.")
			chemid_results = self.handle_no_chemaxon(chemical, request_post)
			if 'error' in chemid_results:
				return chemid_results  # already wrapped for frontend
			if chemid_results.get('smiles'):
				chemical = chemid_results['smiles']
				logging.info("Setting chemical to ACTORWS SMILES to continue chem info routine.")
			elif chemid_results.get('casrn'):
				chemical = chemid_results['casrn']
				logging.info("Setting chemical to ACTORWS casrn to continue chem info routine.")


		# Checks chemical to make sure it's not actually an acronym instead of smiles (e.g., PFAS):
		if chem_type.get('type') == 'smiles' or chem_type.get('type') == 'smarts':
			converted_smiles = self.smiles_name_check(chemical)
			# Switches chem type to "name" if smiles was actually an acronym:
			if converted_smiles:
				# Was actually a name and successfully converted to smiles
				logging.info("Setting chemical to converted_smiles: {}".format(converted_smiles))
				chemical = converted_smiles

		# Uses name form of C, CC, and CCC SMILES:
		if chemical in list(self.carbon_anomolies.keys()):
			chem_type['type'] = "name"
			chemical = self.carbon_anomolies[chemical]
		
		# # ACTORWS requests handling for getting DSSTOX data
		# if chem_type.get('type') == 'CAS#':
		# 	# Gets dsstox results using user-entered CAS#:
		# 	dsstox_results = self.ccte_epa_obj.get_chemical_results(chemical, "casrn")
		# 	_actor_results.update(dsstox_results)
		# elif chem_type.get("type") == "name":
		# 	dsstox_results = self.ccte_epa_obj.get_chemical_results(chemical, "name")
		# 	_actor_results.update(dsstox_results)
		# else:
		# 	pass  # TODO: Use name or cas from ChemAxon to get CCTE data

		# # Returns dsstox substance ID if that's all that's needed,
		# # which is used as the DB key for the chem-info document:
		# if only_dsstox:
		# 	return dsstox_results.get('data', {})

		# # Uses SMILES from actorws if it's there, or user's smiles if not, then jchem smiles if the previous are false:
		# if chem_type['type'] == 'smiles':
		# 	orig_smiles = chemical  # use user-entered smiles as orig_siles
		# elif 'smiles' in _actor_results.get('data', {}) and _actor_results['data']['smiles']:
		# 	orig_smiles = _actor_results['data']['smiles']  # use actorws smiles as orig_smiles
		# 	logging.info("Using actorws smiles as original smiles..")
		# else:
		# 	logging.info("smiles not in user request or actorws results, getting from jchem ws..")
		# 	orig_smiles = self.calc_obj.convertToSMILES({'chemical': chemical}).get('structure')
				# Uses SMILES from actorws if it's there, or user's smiles if not, then jchem smiles if the previous are false:
		if chem_type['type'] == 'smiles':
			orig_smiles = chemical  # use user-entered smiles as orig_siles
		else:
			logging.info("smiles not in user request or actorws results, getting from jchem ws..")
			orig_smiles = self.calc_obj.convertToSMILES({'chemical': chemical}).get('structure')

		# Gets filtered SMILES:
		try:
			filtered_smiles = self.smiles_filter_obj.filterSMILES(orig_smiles, is_node=request_post.get('is_node'))			
			if isinstance(filtered_smiles, dict) and 'error' in filtered_smiles:
				response_obj = {}
				response_obj['status'] = False
				response_obj['request_post'] = request_post
				response_obj['error'] = filtered_smiles['error']
				return response_obj
		except Exception as e:
			logging.warning("Error filtering SMILES: {}".format(e))
			response_obj = {}
			response_obj['status'] = False
			response_obj['error'] = "Cannot process chemical"
			response_obj['request_post'] = request_post
			return response_obj

		# Gets chemical details from jchem ws:
		jchem_response = self.calc_obj.getChemDetails({'chemical': filtered_smiles})
		
		# Creates molecule object with jchem response:
		molecule_obj = Molecule().createMolecule(chemical, orig_smiles, jchem_response, get_sd)

		# Sets 'smiles' (main chemical key for pchem requests, etc.) to CTS standardized smiles:
		molecule_obj['smiles'] = filtered_smiles


		# # ACTORWS requests handling for getting DSSTOX data
		# if chem_type.get('type') == 'CAS#':
		# 	# Gets dsstox results using user-entered CAS#:
		# 	dsstox_results = self.ccte_epa_obj.get_chemical_results(chemical, "casrn")
		# elif chem_type.get("type") == "name":
		# 	dsstox_results = self.ccte_epa_obj.get_chemical_results(chemical, "name")
		# else:
		# 	dsstox_results = self.ccte_epa_obj.get_chemical_results(molecule_obj["preferredName"], "name")
		# _actor_results.update(dsstox_results)

		# Public CCTE requests handling for getting DSSTOX data
		ccte_results = self.ccte_obj.make_search_request(molecule_obj["preferredName"])

		if ccte_results:
			_actor_results.update(ccte_results)
			
		logging.warning("CCTE RESULTS: {}".format(ccte_results))

		# Returns dsstox substance ID if that's all that's needed,
		# which is used as the DB key for the chem-info document:
		if only_dsstox:
			# return dsstox_results.get('data', {})
			return _actor_results.get('data', {})


		cas_list = self.make_cas_request(filtered_smiles)  # gets CAS from cactus.nci.nih.gov (deprecated in jchemws)

		molecule_obj['cas'] = cas_list


		has_carbon = self.smiles_filter_obj.check_for_carbon(filtered_smiles)
		if not has_carbon and is_node:
			molecule_obj['has_carbon'] = False
		else:
			molecule_obj['has_carbon'] = True


		# # Sets 'smiles' (main chemical key for pchem requests, etc.) to CTS standardized smiles:
		# molecule_obj['smiles'] = filtered_smiles

		# Replaces certain keys in molecule_obj with actorws values:
		for key, val in _actor_results.get('data', {}).items():
			if key != 'iupac' and key != 'smiles':
				molecule_obj[key] = val

		# Fills any empty keys with "N/A" for values:
		for key in self.actorws_obj.dsstox_result_keys:
			if key not in molecule_obj:
				molecule_obj.update({key: "N/A"})  # fill in any missed data from actorws with "N/A"

		# Adds popup image with cheminfo table if it's a gentrans product (i.e., node):
		# if is_node or db_handler.is_connected:
		if is_node:
			molecule_obj.update({'node_image': self.calc_obj.nodeWrapper(filtered_smiles, self.calc_obj.tree_image_height, self.calc_obj.tree_image_width, self.calc_obj.image_scale, self.calc_obj.metID,'svg', True)})
			molecule_obj.update({
				'popup_image': self.calc_obj.popupBuilder(
					{"smiles": filtered_smiles}, 
					self.calc_obj.metabolite_keys, 
					"{}".format(request_post.get('id')),
					"Metabolite Information", True)
			})

		wrapped_post = {}
		wrapped_post['status'] = True  # 'metadata': '',
		wrapped_post['data'] = molecule_obj
		wrapped_post['request_post'] = request_post

		return wrapped_post

	def handle_no_chemaxon(self, chemical, request_post):
		"""
		Returns data for ACTORWS only if chemaxon
		isn't available or can't recognize the chemical.
		"""
		molecule_obj = {}
		chemid_results = self.get_chemid_from_actorws(chemical)
		if not chemid_results or not chemid_results.get('smiles'):
			response_obj = {}
			response_obj['status'] = False
			response_obj['request_post'] = request_post
			response_obj['error'] = "Cannot find data for chemical"
			return response_obj
		# remaps keys to cts key names:
		for key, val in chemid_results.items():
			cts_key = self.actorws_obj.chemid_keys_map.get(key)
			if not cts_key:
				continue
			molecule_obj[cts_key] = val
		return molecule_obj

	def smiles_name_check(self, chemical):
		"""
		Known as "the PFOS problem," which is an example chemical of
		an issue where the chemical name is interpretted by JchemWS
		as a SMILES. It tries to convert the chemical into a SMILES, which
		should trigger an error if it actually is one.

		Returns: (True, actual SMILES from JchemWS) if chemical was actual a name,
		(False, original smiles from input) if chemical was actually a smiles.
		"""
		converted_name_response = self.calc_obj.get_smiles_from_name(chemical)
		if converted_name_response.get('smiles') and not 'error' in converted_name_response:
			# if valid, assume chemical was intended to be 'name' instead of 'smiles'..
			# return True
			return converted_name_response['smiles']
		else:
			# if an error was thrown, it was actually smiles so returns original version:
			return None

	def get_chemid_from_actorws(self, chemical, chem_type_name=None):
		_gsid = None
		_smiles_from_mrv = False
		_name_or_smiles = chem_type_name in ['name', 'common', 'smiles', 'systematic']  # bool for chemical in name/common or smiles format
		# If user drew a chemical, get SMILES of chemical from Jchem WS..
		if chem_type_name == 'mrv':
			logging.info("Getting SMILES from jchem web services..")
			response = self.calc_obj.convertToSMILES({'chemical': chemical})
			chemical = response['structure']
			logging.info("SMILES of drawn chemical: {}".format(chemical))
			_smiles_from_mrv = True
		# NOTE: Should be name or smiles, but tries to anyway in case chem type was unknown:
		logging.info("Getting gsid from actorws chemicalIdentifier..")
		chemid_results = self.actorws_obj.get_chemid_results(chemical)  # obj w/ keys calc, prop, data
		return chemid_results

	def make_cas_request(self, smiles):
		"""
		Manually gets CAS list, which used to work with
		Jchem Web Services.
		"""
		try:
			url = self.cas_url.format(requests.utils.quote(smiles))  # encoding smiles for url
			response = requests.get(url, verify=False, timeout=5)
			if response.status_code != 200:
				return "N/A"
			if '<html>' in response.content.decode('utf-8'):
				return "N/A"
			return response.content.decode('utf-8').replace('\n', ', ')  # returns curated CAS list
		except Exception as e:
			logging.warning("Exception making CAS request: {}".format(e))
			return "N/A"

	def check_structure_request(self, chemical):
		"""
		Runs Jchem WS /structureChecker endpoint.
		Initially implemented for validating SMILES
		(e.g., "ccc" should be invalid).
		"""
		url = self.calc_obj.jchem_server_url + self.calc_obj.checker_endpoint
		post = {
			'structure': chemical
		}
		
		results = self.calc_obj.web_call(url, post)

		if not results.get('valid'):
			logging.warning("STRUCTURE CHECK INVALID: {}".format(results))
			return {
				'error': "error requesting structure checker"
			}
		return self.is_valid_aromaticity(results)

	def is_valid_aromaticity(self, check_struct_results):
		"""
		Checks aromaticity results for errors from Jchem WS /structureChecker.
		"""
		assert isinstance(check_struct_results, dict), \
			"is_valid_structure input must be dict of structure results"

		if not "aromaticity" in check_struct_results:
			return {"valid": True}
		elif "description" in check_struct_results["aromaticity"]:
			# returns jchem error to user
			return {"error": check_struct_results["aromaticity"]["description"]}
		else:
			return {"error": "structure is not valid"}

	def check_name_smiles_map(self, chemical):
		"""
		Checks user chemical against map of name-smiles.
		"""
		if chemical.lower() in list(self.chem_names_smiles_map.keys()):
			logging.info("Chemical matches chem_names_smiles_map.")
			return self.chem_names_smiles_map[chemical.lower()]
		else:
			return chemical
