import unittest
import json
import os
import inspect
import datetime
import logging
import sys
from numpy import testing as npt
from tabulate import tabulate
from unittest.mock import Mock, patch
from nose.tools import assert_is_not_none

# local requirements (running pytest at qed level):
_path = os.path.dirname(os.path.abspath(__file__))
if 'cts_celery' in _path:
	from qed.cts_celery.cts_calcs.calculator_chemaxon import JchemCalc
elif 'cts_app' in _path:
	from qed.cts_app.cts_calcs.calculator_chemaxon import JchemCalc

from qed.temp_config.set_environment import DeployEnv



class TestChemaxonCalculator(unittest.TestCase):
	"""
	Unit test class for calculator_chemaxon module, and maybe
	all the calculator_chemaxon depending on how bloated this gets.
	"""

	print("cts calculator_chemaxon unittests conducted at " + str(datetime.datetime.today()))

	def setUp(self):
		"""
		Setup routine for Kabam unit tests.
		:return:
		"""

		# Sets up runtime environment:
		runtime_env = DeployEnv()
		runtime_env.load_deployment_environment()

		# Test chemical if not using mock request json:
		self.test_chemical = "aspirin"
		self.test_smiles = "CC(=O)OC1=C(C=CC=C1)C(O)=O"  # smiles version of aspirin

		# Test inputs for nodeWrapper's smilesToImage function call:
		self.node_wrapper_test_input = {
			'smiles': self.test_smiles,
			'height': None,
			'width': None,
			'scale': None,
			'key': None,
			'img_type': 'svg',
			'isProduct': None
		}

		self.popup_builder_test_input = {
			'smiles': self.test_smiles,
			'iupac': "",
			'formula': "",
			'mass': "",
			'exactMass': "",
			'html': ""
		}

		# Defines filename structure for example JSON results from jchem WS:
		self.filename_structure = "mock_json/calculator_chemaxon_result_{}.json"

		self.calc_obj = JchemCalc()



	def tearDown(self):
		"""
		Teardown routine for Kabam unit tests.
		:return:
		"""
		pass
		# teardown called after each test
		# e.g. maybe write test results to some text file



	def get_example_result_json(self, prop):
		"""
		Gets .json file of example Jchem WS result for unit testing.
		"""
		# filename_structure = "jchem_unittest_object_{}.json"  # {} = property (e.g., pka, tautomer)
		
		filename = self.filename_structure.format(prop)
		
		project_root = os.path.abspath(os.path.dirname(__file__))
		filename = os.path.join(project_root, filename)

		filein = open(filename, 'r')
		file_data = filein.read()
		filein.close()

		return json.loads(file_data)



	@patch('qed.cts_app.cts_calcs.calculator_chemaxon.requests.post')
	@patch('qed.cts_app.cts_calcs.calculator_chemaxon.JchemProperty.getJchemPropData')
	@patch('qed.cts_app.cts_calcs.calculator_chemaxon.SMILESFilter.parseSmilesByCalculator')
	def test_data_request_handler(self, smiles_filter_mock, pchem_mock, speciation_mock):
		"""
		Testing chemaxon calculator module's data_request_handler function.
		"""

		print(">>> Running calculator data_request_handler unit test..")

		smiles_filter_mock.return_value = self.test_smiles  # sets mock val for smiles filter call
		pchem_mock.return_value = {'data': None}  # mock data for pchem request
		speciation_mock.return_value = {'data': None}  # mock data for speciation request

		test_input = {
			'method': None,
			'node': None
		}

		expected_result = self.calc_obj.pchem_request
		response = self.calc_obj.data_request_handler(test_input)

		try:
			self.assertEqual(response['chemical'], expected_result['chemical'])

		finally:
			tab = [response, expected_result]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return