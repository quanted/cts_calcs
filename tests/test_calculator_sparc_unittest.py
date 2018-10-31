import unittest
import json
import os
import inspect
import datetime
import logging
import sys
import requests
from numpy import testing as npt
from tabulate import tabulate
from unittest.mock import Mock, patch
from nose.tools import assert_is_not_none

# local requirements (running pytest at qed level):
_path = os.path.dirname(os.path.abspath(__file__))
if 'cts_celery' in _path:
	from qed.cts_celery.cts_calcs.calculator_sparc import SparcCalc
elif 'cts_app' in _path:
	from qed.cts_app.cts_calcs.calculator_sparc import SparcCalc

from qed.temp_config.set_environment import DeployEnv



class TestSparcCalculator(unittest.TestCase):
	"""
	Unit test class for calculator_sparc module, and maybe
	all the calculator_sparc depending on how bloated this gets.
	"""

	print("cts calculator_sparc unittests conducted at " + str(datetime.datetime.today()))

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
		self.filename_structure = "mock_json/calculator_sparc_result_{}.json"

		self.calc_obj = SparcCalc()



	def tearDown(self):
		"""
		Teardown routine for Kabam unit tests.
		:return:
		"""
		pass
		# teardown called after each test
		# e.g. maybe write test results to some text file



	def test_validate_response(self):
		"""
		Testing Measured's validate_response function.
		"""

		expected_result = True  # expected function result

		response_obj = requests.Response()
		response_obj.status_code = 200

		response = self.calc_obj.validate_response(response_obj)

		try:
			self.assertEqual(response, expected_result)

		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_makeDataRequest(self):
		"""
		Testing Measured's makeDataRequest function.
		"""

		expected_result = {'test': True} # expected function result

		with patch('qed.cts_app.cts_calcs.calculator_epi.EpiCalc.request_logic') as service_mock:
			service_mock.return_value = expected_result
			response = self.calc_obj.makeDataRequest()

		try:
			self.assertDictEqual(response, expected_result)

		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	@patch('qed.cts_app.cts_calcs.calculator_sparc.SparcCalc.validate_response')
	@patch('qed.cts_app.cts_calcs.calculator_sparc.requests.post')	
	def test_request_logic(self, request_mock, validate_mock):
		"""
		Testing Measured's request_logic function.
		"""

		expected_result = {'test': True}  # expected result from test function

		# Mock result for function call in test function:
		request_mock.return_value.content = json.dumps(expected_result)

		validate_mock.return_value = True  # sets mock result for function call in test function

		response = self.calc_obj.request_logic("", {})  # calls test function

		try:
			self.assertDictEqual(response, expected_result)

		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	@patch('qed.cts_app.cts_calcs.calculator_sparc.SparcCalc.makeDataRequest')
	@patch('qed.cts_app.cts_calcs.calculator_sparc.SMILESFilter.parseSmilesByCalculator')
	def test_data_request_handler(self, smiles_filter_mock, request_mock):
		"""
		Testing sparc calculator module's data_request_handler function.
		"""

		print(">>> Running calculator data_request_handler unit test..")

		expected_result = self.calc_obj.pchem_request
		expected_result['chemical'] = self.test_smiles

		smiles_filter_mock.return_value = self.test_smiles  # sets mock val for smiles filter call
		request_mock.return_value = {'test': True}

		test_input = {
			'method': None,
			'prop': 'melting_point',
			'calc': "sparc",
			'chemical': self.test_smiles
		}

		response = self.calc_obj.data_request_handler(test_input)

		try:
			self.assertEqual(response['chemical'], expected_result['chemical'])

		finally:
			tab = [response, expected_result]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return