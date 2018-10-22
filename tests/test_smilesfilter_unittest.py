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
	from qed.cts_celery.cts_calcs.smilesfilter import SMILESFilter
elif 'cts_app' in _path:
	from qed.cts_app.cts_calcs.smilesfilter import SMILESFilter

from qed.temp_config.set_environment import DeployEnv



class TestSmilesFilter(unittest.TestCase):
	"""
	Unit test class for calculators module, and maybe
	all the calculators depending on how bloated this gets.
	"""

	print("cts calculators unittests conducted at " + str(datetime.datetime.today()))

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

		# Defines filename structure for example JSON results from jchem WS:
		self.filename_structure = "mock_json/smilesfilter_result_{}.json"

		self.smilesfilter_obj = SMILESFilter()



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



	def test_is_valid_smiles(self):
		"""
		Testing smilesfilter module is_valid_smiles() function.
		"""

		print(">>> Running smilesfilter is_valid_smiles unit test..")

		mock_json = {"result": "true"}  # expected response from ctsws

		expected_result = True  # expected result from smilesfilter test function

		with patch('qed.cts_app.cts_calcs.smilesfilter.requests.post') as service_mock:

			service_mock.return_value.content = json.dumps(mock_json)  # sets expected result from ctsws request

			response = self.smilesfilter_obj.is_valid_smiles(self.test_smiles)

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_singleFilter(self):
		"""
		Testing smilesfilter module singleFilter() function.
		"""

		print(">>> Running smilesfilter singleFilter unit test..")

		# Test function input dict:
		test_input = {
			'smiles': self.test_smiles,
			'action': "neutralize"
		}

		# Expected response from test function after mock web_call:
		expected_result = {
			"results": [self.test_smiles],
			"status": "success"
		}

		with patch('qed.cts_app.cts_calcs.smilesfilter.Calculator.web_call') as service_mock:

			service_mock.return_value = expected_result  # sets expected result from ctsws request

			response = self.smilesfilter_obj.singleFilter(test_input)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	@patch('qed.cts_app.cts_calcs.smilesfilter.Tautomerization')  # tautomer_mock
	@patch('qed.cts_app.cts_calcs.smilesfilter.Calculator.web_call')  # ctsws_call_2
	@patch('qed.cts_app.cts_calcs.smilesfilter.Calculator.web_call')  # ctsws_call_1
	@patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.is_valid_smiles')  # validity check mock
	def test_filterSMILES(self, validity_check_mock, ctsws_call_1, ctsws_call_2, tautomer_mock):
		"""
		Testing smilesfilter module filterSMILES() function.

		NOTE: This service makes 3 requests to servers, so there are
		3 patches that mock those responses. Order of patches and input
		to this test function are ordered bottom to top.
		"""

		print(">>> Running smilesfilter filterSMILES unit test..")

		# Expected value from first web_call mock:
		expected_mock_val_1 = {
		    "results": [
		        self.test_smiles,
		        self.test_smiles
		    ],
		    "status": "success"
		}

		# Expected value from 2nd web_call mock:
		expected_mock_val_2 = {
		    "results": [
		        self.test_smiles
		    ],
		    "status": "success"
		}

		# Setting up tautomer mock:
		self.filename_structure = "mock_json/smilesfilter_result_{}.json"  # setting filename structure to jchem_properties json results
		expected_tautomer_mock_val = self.get_example_result_json("tautomerization")  # opens tautomer result json for mocking request

		# Test function input dict:
		test_input = {
			'smiles': self.test_smiles,
			'action': "neutralize"
		}

		expected_result = self.test_smiles  # expecting same smiles after filter for this example

		# Sets mock return values for function calls within test function:00
		validity_check_mock.return_value = True
		ctsws_call_1.return_value = expected_mock_val_1  # sets mock val for 1st ctsws smiles filter request
		ctsws_call_2.return_value = expected_mock_val_2  # sets mock val for 2nd ctsws smiles filter request
		tautomer_mock.return_value.make_data_request.return_value = expected_tautomer_mock_val  # sets mock val for jchemws taut request
		tautomer_mock.return_value.results = expected_tautomer_mock_val

		response = self.smilesfilter_obj.filterSMILES(self.test_smiles)  # runs test function

		try:
			# Compares function response with expected:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_checkMass(self):
		"""
		Testing smilesfilter module checkMass() function.
		"""

		print(">>> Running smilesfilter checkMass unit test..")

		# Setting up mass mock:
		self.filename_structure = "mock_json/calculator_result_{}.json"  # setting filename structure to jchem_properties json results
		expected_mass_mock_val = self.get_example_result_json("getmass")  # opens mass result json for mocking request

		# Expected response from test function:
		expected_result = True

		with patch('qed.cts_app.cts_calcs.smilesfilter.Calculator.getMass') as service_mock:

			service_mock.return_value = expected_mass_mock_val  # sets expected result from getMass request

			response = self.smilesfilter_obj.checkMass(self.test_chemical)

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_clearStereos(self):
		"""
		Testing smilesfilter module clearStereos() function.
		"""

		print(">>> Running smilesfilter clearStereos unit test..")

		# Expected result from mock function called by test function:
		mock_json = {
			"results": [self.test_smiles],
			"status": "success"
		}

		# Expected response from test function:
		expected_result = self.test_smiles

		with patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.singleFilter') as service_mock:

			service_mock.return_value = mock_json  # sets expected result from getMass request

			response = self.smilesfilter_obj.clearStereos(self.test_smiles)
			response = response[0]  # expecting one smiles in list after test function call

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_transformSMILES(self):
		"""
		Testing smilesfilter module transformSMILES() function.
		"""

		print(">>> Running smilesfilter transformSMILES unit test..")

		# Expected result from mock function called by test function:
		mock_json = {
			"results": [self.test_smiles],
			"status": "success"
		}

		# Expected response from test function:
		expected_result = self.test_smiles

		with patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.singleFilter') as service_mock:

			service_mock.return_value = mock_json  # sets expected result from getMass request

			response = self.smilesfilter_obj.transformSMILES(self.test_smiles)
			response = response[0]  # expecting one smiles in list after test function call

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_untransformSMILES(self):
		"""
		Testing smilesfilter module untransformSMILES() function.
		"""

		print(">>> Running smilesfilter untransformSMILES unit test..")

		# Expected result from mock function called by test function:
		mock_json = {
			"results": [self.test_smiles],
			"status": "success"
		}

		# Expected response from test function:
		expected_result = self.test_smiles

		with patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.singleFilter') as service_mock:

			service_mock.return_value = mock_json  # sets expected result from getMass request

			response = self.smilesfilter_obj.untransformSMILES(self.test_smiles)
			response = response[0]  # expecting one smiles in list after test function call

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	@patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.untransformSMILES')
	@patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.clearStereos')
	@patch('qed.cts_app.cts_calcs.smilesfilter.SMILESFilter.checkMass')
	def test_parseSmilesByCalculator(self, mass_mock, stereos_mock, untransform_mock):
		"""
		Testing smilesfilter module parseSmilesByCalculator() function.
		"""

		print(">>> Running smilesfilter parseSmilesByCalculator unit test..")

		test_input = "epi"  # calc input for test function

		# Setting mock return vals for test function:
		mass_mock.return_value = True
		stereos_mock.return_value = [self.test_smiles]
		untransform_mock.return_value = [self.test_smiles]

		# Expected response from test function:
		expected_result = self.test_smiles

		response = self.smilesfilter_obj.parseSmilesByCalculator(self.test_smiles, test_input)	

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return