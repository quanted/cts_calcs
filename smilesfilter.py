import requests
import json
import logging
from .calculator import Calculator
from .jchem_properties import Tautomerization



class SMILESFilter(object):
	"""
	This is the smilesfilter.py module as a class and
	clumped in with other classes related to chem info.
	"""

	def __init__(self):
		self.max_weight = 1500  # max weight [g/mol] for epi, test, and sparc
		self.excludestring = {".","[Ag]","[Al]","[As","[As+","[Au]","[B]","[B-]","[Br-]","[Ca]",
						"[Ca+","[Cl-]","[Co]","[Co+","[Fe]","[Fe+","[Hg]","[K]","[K+","[Li]",
						"[Li+","[Mg]","[Mg+","[Na]","[Na+","[Pb]","[Pb2+]","[Pb+","[Pt]",
						"[Sc]","[Si]","[Si+","[SiH]","[Sn]","[W]"}
		self.return_val = {
			"valid" : False,
			"smiles": "",
			"processedsmiles" : ""
		}

	def is_valid_smiles(self, smiles):

		if any(x in smiles for x in self.excludestring):
			return False

		return True


	def singleFilter(self, request_obj):
		"""
		Calls single EFS Standardizer filter
		for filtering SMILES
		"""
		smiles = request_obj.get('smiles')
		action = request_obj.get('action')
		post_data = {
			"structure": smiles,
			"actions": [
				action
			]
		}
		calc = Calculator()
		url = calc.efs_server_url + calc.efs_standardizer_endpoint
		return calc.web_call(url, post_data)


	def filterSMILES(self, smiles):
		"""
		cts ws call to jchem to perform various
		smiles processing before being sent to
		p-chem calculators
		"""
		calc_object = Calculator()

		logging.info("{} is being processed by cts SMILES filter..".format(smiles))
		logging.info("Checking chemical for any metals..")

		if not self.is_valid_smiles(smiles):
			logging.warning("User chemical contains metals, sending error to client..")
			# raise Exception({'data': "Chemical cannot contain metals.."})
			raise ValueError("Chemical cannot contain metals..")

		# Updated approach (todo: more efficient to have CTSWS use major taut instead of canonical)
		# 1. CTSWS actions "removeExplicitH" and "transform".
		url = calc_object.efs_server_url + calc_object.efs_standardizer_endpoint
		post_data = {
			'structure': smiles,
			'actions': [
				"removeExplicitH",
				"transform"
			]
		}
		response = calc_object.web_call(url, post_data)

		logging.info("1. Removing explicit H, then transforming..")
		logging.info("request to jchem: {}".format(post_data))
		logging.info("response from jchem: {}".format(response))

		filtered_smiles = response['results'][-1] # picks last item, format: [filter1 smiles, filter1 + filter2 smiles]

		logging.info("filtered smiles so far: {}".format(filtered_smiles))
		
		# 2. Get major tautomer from jchem:
		taut_obj = Tautomerization()
		taut_obj.postData.update({'calculationType': 'MAJOR'})
		taut_obj.make_data_request(filtered_smiles, taut_obj)

		logging.info("2. Obtaining major tautomer from {}".format(filtered_smiles))
		logging.info("request to jchem: {}".format(taut_obj.postData))
		logging.info("response from jchem: {}".format(taut_obj.results))

		# todo: verify this is major taut result smiles, not original smiles for major taut request...
		major_taut_smiles = None
		try:
			major_taut_smiles = taut_obj.results['result']['structureData']['structure']
		except KeyError as e:
			logging.info("Jchem error requesting major tautomer from {}..".format(filtered_smiles))
			logging.info("Using smiles {} for next step..".format(filtered_smiles))

		if major_taut_smiles:
			logging.info("Major tautomer found: {}.. Using as filtered smiles..".format(major_taut_smiles))
			filtered_smiles = major_taut_smiles

		# 3. Using major taut smiles for final "neutralize" filter:
		post_data = {
			'structure': filtered_smiles, 
			'actions': [
				"neutralize"
			]
		}
		response = calc_object.web_call(url, post_data)

		logging.info("3. Neutralizing smiles {}".format(filtered_smiles))
		logging.info("request to jchem: {}".format(post_data))
		logging.info("response from jchem: {}".format(response))

		final_smiles = response['results'][-1]
		logging.info("smiles results after cts filtering: {}".format(response.get('results')))
		logging.info("FINAL FITERED SMILES: {}".format(final_smiles))

		return final_smiles


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
				# raise "Structure too large, must be < 1500 g/mol.."
				raise Exception({'data': "structure too large"})

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
				raise {'data': "error filtering chemical"}

		# 4. Check for metals and stuff (square brackets):
		if calculator == 'epi' or calculator == 'measured':
			if '[' in filtered_smiles or ']' in filtered_smiles:
				# bubble up to calc for handling error
				# raise Exception("{} cannot process metals..".format(calculator))
				raise Exception({'data': "cannot process metals or charges"})

		return filtered_smiles