__author__ = "np"

import requests
import logging
import json
from .calculator_metabolizer import MetabolizerCalc
from .actorws import ACTORWS
from .smilesfilter import SMILESFilter




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
		self.cas = ''
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
					if key == 'structureData' and get_structure_data == None:
						pass
					elif key == 'cas':
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
		self.wrapped_post = {
			'status': False,
			'data': None,
			'request_post': None
		}


	def get_cheminfo(self, request_post):
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

		actorws = ACTORWS()
		smiles_filter = SMILESFilter()
		calc = MetabolizerCalc()  # note: inherits Calculator class as well
		_actor_results = {}
		orig_smiles = ""  # initial SMILES pre CTS filter

		# Determine chemical type from user (e.g., smiles, cas, name, etc.):
		chem_type = calc.get_chemical_type(chemical)

		if 'error' in chem_type:
			response_obj = self.wrapped_post
			response_obj['error'] = chem_type['error']
			response_obj['request_post'] = request_post
			return response_obj

		logging.info("Incoming chemical to CTS standardizer: {}".format(chemical))
		logging.info("Chemical type: {}".format(chem_type))

		if chem_type.get('type') == 'smiles':
			# try conversion to smiles assuming chemical is name to break any ties if chem valid name and smiles
			# inspired by the PFOS problem..
			converted_name_response = calc.get_smiles_from_name(chemical)
			if converted_name_response.get('smiles') and not 'error' in converted_name_response:
				# if valid, assume chemical was intended to be 'name' instead of 'smiles'..
				orig_smiles = converted_name_response.get('smiles')  # used converted smiles from name

		_gsid = self.get_chemid_from_actorws(chemical, chem_type.get('type'), actorws)

		# Should be CAS# or have gsid from chemid by this point..
		if _gsid or chem_type['type'] == 'CAS#':
			id_type = 'CAS#'
			if _gsid:
				chem_id = _gsid  # use gsid for ACTORWS request
				id_type = 'gsid'
			else:
				chem_id = chemical  # use CAS# for ACTORWS request
			logging.info("Getting results from actorws dsstox..")
			dsstox_results = actorws.get_dsstox_results(chem_id, id_type)  # keys: smiles, iupac, preferredName, dsstoxSubstanceId, casrn 
			_actor_results.update(dsstox_results)

		# If user enters something other than SMILES, use actorws smiles for orig_smiles
		if orig_smiles:
			pass
		elif chem_type['type'] == 'smiles':
			orig_smiles = chemical  # use user-entered smiles as orig_siles
		elif 'smiles' in _actor_results.get('data', {}):
			orig_smiles = _actor_results['data']['smiles']  # use actorws smiles as orig_smiles
			logging.info("Using actorws smiles as original smiles..")
		else:
			logging.info("smiles not in user request or actorws results, getting from jchem ws..")
			orig_smiles = calc.convertToSMILES({'chemical': chemical}).get('structure')

		try:
			filtered_smiles = smiles_filter.filterSMILES(orig_smiles)
		except ValueError as e:
			logging.warning("Error filtering SMILES: {}".format(e))
			response_obj = self.wrapped_post
			response_obj['error'] = "Chemical cannot contain metals.."
			response_obj['request_post'] = request_post
			return response_obj
		except Exception as e:
			logging.warning("Error filtering SMILES: {}".format(e))
			response_obj = self.wrapped_post
			response_obj['error'] = "Cannot process chemical.."
			response_obj['request_post'] = request_post
			return response_obj


		logging.info("Original smiles before cts filtering: {}".format(orig_smiles))
		logging.info("Filtered SMILES: {}".format(filtered_smiles))

		jchem_response = calc.getChemDetails({'chemical': filtered_smiles})  # get chemical details

		molecule_obj = Molecule().createMolecule(chemical, orig_smiles, jchem_response, get_sd)

		# Loop _actor_results, replace certain keys in molecule_obj with actorws vals:
		for key, val in _actor_results.get('data', {}).items():
			if key != 'iupac' and key != 'smiles':
				# using chemaxon 'iupac' instead of actorws
				molecule_obj[key] = val  # replace or add any values from chemaxon deat

		for key in actorws.dsstox_result_keys:
			if key not in molecule_obj and key != 'preferredName':
				# Note: using preferredName from chemaxon instead..
				molecule_obj.update({key: "N/A"})  # fill in any missed data from actorws with "N/A"

		if is_node:
			molecule_obj.update({'node_image': calc.nodeWrapper(filtered_smiles, calc.tree_image_height, calc.tree_image_width, calc.image_scale, calc.metID,'svg', True)})
			molecule_obj.update({
				'popup_image': Calculator().popupBuilder(
					{"smiles": filtered_smiles}, 
					calc.metabolite_keys, 
					"{}".format(request_post.get('id')),
					"Metabolite Information")
			})

		wrapped_post = self.wrapped_post
		wrapped_post['status'] = True  # 'metadata': '',
		wrapped_post['data'] = molecule_obj
		wrapped_post['request_post'] = request_post
		return wrapped_post


	def get_chemid_from_actorws(self, chemical, chem_type_name, actorws):
		_gsid = None
		_smiles_from_mrv = False
		_name_or_smiles = chem_type_name in ['name', 'common', 'smiles']  # bool for chemical in name/common or smiles format

		# If user drew a chemical, get SMILES of chemical from Jchem WS..
		if chem_type_name == 'mrv':
			logging.info("Getting SMILES from jchem web services..")
			response = calc.convertToSMILES({'chemical': chemical})
			chemical = response['structure']
			logging.info("SMILES of drawn chemical: {}".format(chemical))
			_smiles_from_mrv = True

		if _name_or_smiles or _smiles_from_mrv:
			logging.info("Getting gsid from actorws chemicalIdentifier..")
			chemid_results = actorws.get_chemid_results(chemical)  # obj w/ keys calc, prop, data
			_gsid = chemid_results.get('data', {}).get('gsid')
			logging.info("gsid from actorws chemid: {}".format(_gsid))

		return _gsid