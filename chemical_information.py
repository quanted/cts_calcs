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
		Main function for Chemical Editor. Return chemical info from CCTE/Comptox and
		uses rdkit as well.
		"""

		chemical = request_post.get('chemical')
		get_sd = request_post.get('get_structure_data')  # bool for getting <cml> format image for marvin sketch
		is_node = request_post.get('is_node')  # bool for tree node or not
		orig_smiles = None  # initial SMILES pre CTS filter

		# Checks chemical against chem_name_smiles_map:
		chemical = self.check_name_smiles_map(chemical)

		chem_type = self.determine_chem_type(chemical)

		# Gets initial chem info data from CCTE as well as sets orig_smiles:
		ccte_results = self.get_initial_ccte_data(chemical, chem_type, request_post)
		if "error" in ccte_results:
			return ccte_results

		# Uses name form of C, CC, and CCC SMILES:
		if chemical in list(self.carbon_anomolies.keys()):
			chemical = self.carbon_anomolies[chemical]

		logging.info("Original SMILES: {}".format(ccte_results["data"]["orig_smiles"]))

		# Gets filtered SMILES:
		filtered_smiles = self.make_smilesfilter_request(ccte_results["data"]["orig_smiles"], request_post)
		if "error" in filtered_smiles:
			return filtered_smiles

		logging.info("Filtered SMILES: {}".format(filtered_smiles))

		# Gets detailed chem info from CCTE using dtxsid:
		ccte_detail_results = self.ccte_obj.make_details_request(ccte_results["data"]["dtxsid"])

		final_results = {}
		for key, val in ccte_detail_results.items():
			# remaps keys to cts key names
			cts_key = self.ccte_obj.chemid_keys_map.get(key)
			if not cts_key:
				continue
			final_results[cts_key] = val

		# Returns dsstox substance ID if that's all that's needed,
		# which is used as the DB key for the chem-info document:
		if only_dsstox:
			return ccte_detail_results.get('data', {})

		cas_list = self.make_cas_request(filtered_smiles)  # gets CAS from cactus.nci.nih.gov (deprecated in jchemws)

		final_results['smiles'] = filtered_smiles
		final_results['cas'] = cas_list
		final_results['chemical'] = chemical
		final_results['orig_smiles'] = ccte_results["data"]["orig_smiles"]

		has_carbon = self.smiles_filter_obj.check_for_carbon(filtered_smiles)
		if not has_carbon and is_node:
			final_results['has_carbon'] = False
		else:
			final_results['has_carbon'] = True

		# Adds popup image with cheminfo table if it's a gentrans product (i.e., node):
		if is_node:
			final_results.update({'node_image': self.calc_obj.nodeWrapper(filtered_smiles, self.calc_obj.tree_image_height, self.calc_obj.tree_image_width, self.calc_obj.image_scale, self.calc_obj.metID,'svg', True)})
			final_results.update({
				'popup_image': self.calc_obj.popupBuilder(
					{"smiles": filtered_smiles}, 
					self.calc_obj.metabolite_keys, 
					"{}".format(request_post.get('id')),
					"Metabolite Information", True)
			})

		wrapped_post = {}
		wrapped_post['status'] = True  # 'metadata': '',
		wrapped_post['data'] = final_results
		wrapped_post['request_post'] = request_post

		return wrapped_post

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

	def make_smilesfilter_request(self, orig_smiles, request_post):
		try:
			filtered_smiles = self.smiles_filter_obj.filterSMILES(orig_smiles, is_node=request_post.get('is_node'))			
			if isinstance(filtered_smiles, dict) and 'error' in filtered_smiles:
				response_obj = {}
				response_obj['status'] = False
				response_obj['request_post'] = request_post
				response_obj['error'] = filtered_smiles['error']
				return response_obj
			else:
				return filtered_smiles
		except Exception as e:
			logging.warning("Error filtering SMILES: {}".format(e))
			response_obj = {}
			response_obj['status'] = False
			response_obj['error'] = "Cannot process chemical"
			response_obj['request_post'] = request_post
			return response_obj

	def get_initial_ccte_data(self, chemical, chem_type, request_post):
		if chem_type in ["cas", "dtxsid", "name"]:
			ccte_results = self.ccte_obj.make_search_request(chemical)
			logging.info("ccte_results from make_search_request: {}".format(ccte_results))
			ccte_results["data"]["orig_smiles"] = ccte_results["data"]["smiles"]
			return ccte_results
		elif chem_type == "smiles":
			ccte_results = self.ccte_obj.make_smiles_request(chemical)
			logging.info("ccte_results from make_smiles_request: {}".format(ccte_results))
			ccte_results["data"]["smiles"] = chemical
			ccte_results["data"]["orig_smiles"] = chemical
			return ccte_results
		else:
			logging.error("Chemical type not recognized: {}".format(e))
			response_obj = {}
			response_obj['status'] = False
			response_obj['error'] = "Cannot process chemical"
			response_obj['request_post'] = request_post
			return response_obj

	def determine_chem_type(self, chemical):
		if self.is_cas(chemical):
			logging.info("chem_type: cas")
			return "cas"
		if self.is_dtxsid(chemical):
			logging.info("chem_type: dtxsid")
			return "dtxsid"
		if self.is_valid_smiles(chemical):
			logging.info("chem_type: smiles")
			return "smiles"
		logging.info("chem_type: name")
		return "name"

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


