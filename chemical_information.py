__author__ = "np"

import requests
import logging
import json
import re

from rdkit import Chem
from .calculator_metabolizer import MetabolizerCalc
from .actorws import ACTORWS, CCTE_EPA
from .smilesfilter import SMILESFilter
from .ccte import CCTE



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
			"pfas": "Perfluorooctanoic acid",
			"pfos": "Perfluorooctanesulfonic acid"
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
			{'Preferred Name': workflow_obj.preferredName},
			{'IUPAC': workflow_obj.name},
			{'Formula': workflow_obj.formula},
			{'Preferred CAS': workflow_obj.casrn},
			{'CAS #': workflow_obj.cas},
			{'DTXSID': workflow_obj.dtxsid},
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

		chemical = request_post.get('chemical')
		get_sd = request_post.get('get_structure_data')  # bool for getting <cml> format image for marvin sketch
		is_node = request_post.get('is_node')  # bool for tree node or not
		_actor_results = {}  # dict for actorws results
		_gsid = None
		orig_smiles = None  # initial SMILES pre CTS filter
		is_name = False  # bool for whether smiles was actually acronym

		# # Checks chemical against chem_name_smiles_map:
		chemical = self.check_name_smiles_map(chemical)

		is_cas = self.is_cas(chemical)
		logging.info("is_cas: {}".format(is_cas))

		is_dtxsid = self.is_dtxsid(chemical)
		logging.info("is_dtxsid: {}".format(is_dtxsid))

		is_smiles = self.is_valid_smiles(chemical)
		logging.info("is_smiles: {}".format(is_smiles))

		is_name = False
		if not (is_cas or is_dtxsid or is_smiles):
			is_name = True


		ccte_results = None

		# If not cas, dtxsid, or smiles, assuming name and will try to get SMILES (error assumes invalid structure)
		if is_cas or is_dtxsid or is_name:
			logging.warning("Calling make_search_request")
			# ccte_results = self.ccte_obj.make_search_request(molecule_obj["preferredName"])
			ccte_results = self.ccte_obj.make_search_request(chemical)
			logging.warning("ccte_results: {}".format(ccte_results))

			orig_smiles = ccte_results["data"]["smiles"]

		elif is_smiles:
			logging.warning("Calling make_smiles_request")
			ccte_results = self.ccte_obj.make_smiles_request(chemical)
			ccte_results["data"]["smiles"] = chemical
			orig_smiles = chemical

			logging.warning("ccte_results: {}".format(ccte_results))

		else:
			logging.error("Chemical type not recognized: {}".format(e))
			response_obj = {}
			response_obj['status'] = False
			response_obj['error'] = "Cannot process chemical"
			response_obj['request_post'] = request_post
			return response_obj


		# Uses name form of C, CC, and CCC SMILES:
		if chemical in list(self.carbon_anomolies.keys()):
			chemical = self.carbon_anomolies[chemical]


		logging.info("Original SMILES: {}".format(orig_smiles))


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

		logging.info("Filtered SMILES: {}".format(filtered_smiles))

		ccte_detail_results = self.ccte_obj.make_details_request(ccte_results["data"]["dtxsid"])

		ccte_detail_results_curated = {}

		# remaps keys to cts key names:
		for key, val in ccte_detail_results.items():
			cts_key = self.ccte_obj.chemid_keys_map.get(key)
			if not cts_key:
				continue
			ccte_detail_results_curated[cts_key] = val

		# Returns dsstox substance ID if that's all that's needed,
		# which is used as the DB key for the chem-info document:
		if only_dsstox:
			# return dsstox_results.get('data', {})
			return ccte_detail_results.get('data', {})

		cas_list = self.make_cas_request(filtered_smiles)  # gets CAS from cactus.nci.nih.gov (deprecated in jchemws)

		ccte_detail_results_curated['smiles'] = filtered_smiles
		ccte_detail_results_curated['cas'] = cas_list
		ccte_detail_results_curated['chemical'] = chemical
		ccte_detail_results_curated['orig_smiles'] = orig_smiles


		has_carbon = self.smiles_filter_obj.check_for_carbon(filtered_smiles)
		if not has_carbon and is_node:
			ccte_detail_results_curated['has_carbon'] = False
		else:
			ccte_detail_results_curated['has_carbon'] = True

		# Adds popup image with cheminfo table if it's a gentrans product (i.e., node):
		if is_node:
			ccte_detail_results_curated.update({'node_image': self.calc_obj.nodeWrapper(filtered_smiles, self.calc_obj.tree_image_height, self.calc_obj.tree_image_width, self.calc_obj.image_scale, self.calc_obj.metID,'svg', True)})
			ccte_detail_results_curated.update({
				'popup_image': self.calc_obj.popupBuilder(
					{"smiles": filtered_smiles}, 
					self.calc_obj.metabolite_keys, 
					"{}".format(request_post.get('id')),
					"Metabolite Information", True)
			})

		wrapped_post = {}
		wrapped_post['status'] = True  # 'metadata': '',
		wrapped_post['data'] = ccte_detail_results_curated
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

		logging.info("smiles_name_check converted_name_response: {}".format(converted_name_response))

		if converted_name_response.get('smiles') and not 'error' in converted_name_response:
			# if valid, assume chemical was intended to be 'name' instead of 'smiles'..
			# return True
			logging.info("Received valid response, assuming chemical was intended to be a name instead of smiles.")
			return converted_name_response['smiles']
		else:
			# if an error was thrown, it was actually smiles so returns original version:
			# return None
			logging.info("Assuming chemical was actually smiles and not intended to be a name. Returning original chemical.")
			return chemical

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

	def check_name_smiles_map(self, chemical):
		"""
		Checks user chemical against map of name-smiles.
		"""
		if chemical.lower() in list(self.chem_names_smiles_map.keys()):
			logging.info("Chemical matches chem_names_smiles_map.")
			return self.chem_names_smiles_map[chemical.lower()]
		else:
			return chemical

	def is_valid_smiles(self, chemical):
		try:
			mol = Chem.MolFromSmiles(chemical)
			if mol is None:
				logging.warning("Not valid SMILES.")
				return False
			else:
				return True
		except Exception as e:
			logging.warning("is_valid_smiles exception: {}.\nNot a valid SMILES.".format(e))
			return False
		

	def is_cas(self, chemical):
	    chemical = chemical.strip()
	    if not re.fullmatch(r"\d{2,7}-\d{2}-\d", chemical):
	        return False
	    digits = chemical.replace("-", "")
	    body = digits[:-1]
	    check_digit = int(digits[-1])
	    total = 0
	    for position, digit in enumerate(reversed(body), start=1):
	        total += position * int(digit)
	    return total % 10 == check_digit

	def is_dtxsid(self, chemical):
		chemical = chemical.strip()
		if re.fullmatch(r"DTXSID\d+", chemical, re.IGNORECASE):
			return True
		else:
			return False


