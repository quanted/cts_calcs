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
	from qed.cts_celery.cts_calcs.calculator import Calculator
elif 'cts_app' in _path:
	from qed.cts_app.cts_calcs.calculator import Calculator

from qed.temp_config.set_environment import DeployEnv



class TestCalculators(unittest.TestCase):
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

		# Defines filename structure for example JSON results from jchem WS:
		self.filename_structure = "mock_json/calculator_result_{}.json"

		self.calc_obj = Calculator()



	def tearDown(self):
		"""
		Teardown routine for Kabam unit tests.
		:return:
		"""
		pass
		# teardown called after each test
		# e.g. maybe write test results to some text file

	

	# def create_jchem_object(self, prop):
	# 	pass



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



	def test_gen_jid(self):
		"""
		Testing calculator module's gen_jid function.
		"""

		print(">>> Running calculator gen_jid test..")

		expected_time_format = '%Y%m%d%H%M%S%f'

		expected_results = [20]  # expected length of jid

		# test_json = self.get_example_result_json('genjid')

		jid = self.calc_obj.gen_jid()  # generates a timestamp

		jid_datetime = datetime.datetime.strptime(jid, expected_time_format)  # converts jid to datetime object

		results = [len(jid)]

		try:
			# Comparing expected and results for jid:
			npt.assert_allclose(results, expected_results, rtol=1, atol=0, err_msg='', verbose=True)

			# Checking that jid is valid datetime:
			self.assertIsInstance(jid_datetime, datetime.datetime)

		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	# def test_get_melting_point(self):
	# 	"""
	# 	Testing calculator module's get_melting_point function.
	# 	"""

	# 	print(">>> Running calculator get_melting_point test..")



	def test_get_chem_details(self):
		"""
		Testing calculator module's getChemDetails function, which
		calls JchemWS for chemical information/details.
		"""

		print(">>> Running calculator getChemDetails test..")

		test_chemical = "CC(=O)OC1=C(C=CC=C1)C(O)=O"  # chemical used in expected request and results

		# Gets expected response JSON:
		expected_json = self.get_example_result_json("chem_details")

		# Testing getChemDetails with mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:
			
			# Sets web_call mock to return expected json when getChemDetails calls it:
			service_mock.return_value = expected_json
			
			# Calls getChemDetails, which will use the mock web_call for unit testing:
			response = self.calc_obj.getChemDetails({'chemical': test_chemical})

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_json)
		finally:
			tab = [[response], [expected_json]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_smiles_to_image(self):
		"""
		Testing calculator module's smilesToImage function, which calls
		JchemWS for converting a SMILES to an image.
		"""

		print(">>> Running calculator smilesToImage test..")

		test_chemical = "CC(=O)OC1=C(C=CC=C1)C(O)=O"

		expected_json = self.get_example_result_json("smilestoimage")  # gets expected response json

		# Testing function with a mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:

			# Sets web_call mock to return expected json when test function calls it:
			service_mock.return_value = expected_json

			# Calls test function, which will use the mock web_call for unit testing:
			response = self.calc_obj.smilesToImage({'smiles': test_chemical})

			try:
				# Compares function response with expected json:
				self.assertDictEqual(response, expected_json)
			finally:
				tab = [[response], [expected_json]]
				print("\n")
				print(inspect.currentframe().f_code.co_name)
				print(tabulate(tab, headers='keys', tablefmt='rst'))

