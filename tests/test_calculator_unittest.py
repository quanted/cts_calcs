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



	def get_example_result_html(self, prop):
		"""
		Gets .html file of example Jchem WS result for unit testing.
		"""
		filename_structure = "mock_html/calculator_result_{}.html"
		
		filename = filename_structure.format(prop)
		
		project_root = os.path.abspath(os.path.dirname(__file__))
		filename = os.path.join(project_root, filename)

		filein = open(filename, 'r')
		file_data = filein.read()
		filein.close()

		return file_data



	def test_gen_jid(self):
		"""
		Testing calculator module's gen_jid function.
		"""

		print(">>> Running calculator gen_jid unit test..")

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

	# 	print(">>> Running calculator get_melting_point unit test..")



	def test_getChemDetails(self):
		"""
		Testing calculator module's getChemDetails function, which
		calls JchemWS for chemical information/details.
		"""

		print(">>> Running calculator getChemDetails unit test..")

	# chemical used in expected request and results

		# Gets expected response JSON:
		expected_json = self.get_example_result_json("chem_details")

		# Testing getChemDetails with mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:
			
			# Sets web_call mock to return expected json when getChemDetails calls it:
			service_mock.return_value = expected_json
			
			# Calls getChemDetails, which will use the mock web_call for unit testing:
			response = self.calc_obj.getChemDetails({'chemical': self.test_smiles})

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_json)
		finally:
			tab = [[response], [expected_json]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_smilesToImage(self):
		"""
		Testing calculator module's smilesToImage function, which calls
		JchemWS for converting a SMILES to an image.
		"""

		print(">>> Running calculator smilesToImage unit test..")


		expected_json = self.get_example_result_json("smiles_to_image")  # gets expected response json

		# Testing function with a mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:

			# Sets web_call mock to return expected json when test function calls it:
			service_mock.return_value = expected_json

			# Calls test function, which will use the mock web_call for unit testing:
			response = self.calc_obj.smilesToImage({'smiles': self.test_smiles})

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_json)
		finally:
			tab = [[response], [expected_json]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_convertToSMILES(self):
		"""
		Testing calculator module's convertToSmiles function, which calls
		JchemWS for converting a chemical to a SMILES.
		"""

		print(">>> Running calculator convertToSmiles unit test..")

		expected_json = self.get_example_result_json("convert_to_smiles")  # gets expected response json

		# Testing function with a mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:

			# Sets web_call mock to return expected json when test function calls it:
			service_mock.return_value = expected_json

			# Calls test function, which will use the mock web_call for unit testing:
			response = self.calc_obj.convertToSMILES({'structure': self.test_chemical})

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_json)
		finally:
			tab = [[response], [expected_json]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_getStructInfo(self):
		"""
		Testing calculator module's getStructInfo function, which calls
		calculator's getChemDetails and returns a wrapped dictionary of
		chemical info.
		"""

		print(">>> Running calculator getStructInfo unit test..")

		expected_keys = ['formula', 'iupac', 'mass', 'smiles','exactMass']

		# Mock JSON for getChemDetails, which is called by getStructInfo:
		mock_json = self.get_example_result_json("chem_details")

		# Expected dictionary from getStructInfo:
		expected_dict = {
			'formula': "C9H8O4",
			'iupac': "2-(acetyloxy)benzoic acid",
			'mass': 180.159,
			'exactMass': 180.042258738,
			'smiles': self.test_smiles
		}


		# Testing function with a mock of getChemDetails:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.getChemDetails') as service_mock:

			# Sets getChemDetails mock to return expected json when test function calls it:
			service_mock.return_value = mock_json

			# Calls test function, which will use the getChemDetails mock for unit testing:
			response = self.calc_obj.getStructInfo(self.test_smiles)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_dict)
		finally:
			tab = [[response], [expected_dict]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_getMass(self):
		"""
		Testing calculator module's getMass function, which calls
		JchemWS for a chemical's mass.
		"""

		print(">>> Running calculator getMass unit test..")

		expected_json = self.get_example_result_json("getmass")  # gets expected response json

		# Testing function with a mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:

			# Sets web_call mock to return expected json when test function calls it:
			service_mock.return_value = expected_json

			# Calls test function, which will use the mock web_call for unit testing:
			response = self.calc_obj.getMass({'chemical': self.test_smiles})

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_json)
		finally:
			tab = [[response], [expected_json]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_get_chemical_type(self):
		"""
		Testing calculator module's get_chemical_type function, which calls
		JchemWS for a chemical's format type.
		"""

		print(">>> Running calculator get_chemical_type unit test..")

		expected_json = self.get_example_result_json("get_chemical_type")  # gets expected response json for mock requests.post
		expected_response = {'type': "smiles"}  # expected response from test function

		# Testing function with a mock of requests.post, which is used by test function:
		with patch('qed.cts_app.cts_calcs.calculator.requests.post') as service_mock:

			# Sets requests.post mock to return expected response content when test function calls it:
			service_mock.return_value.content = json.dumps(expected_json)

			# Calls test function, which will use the mock web_call for unit testing:
			response = self.calc_obj.get_chemical_type(self.test_smiles)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_response)
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_get_smiles_from_name(self):
		"""
		Testing calculator module's get_smiles_from_name function, which calls
		JchemWS to get a SMILES from a chemical name.
		"""

		print(">>> Running calculator get_smiles_from_name unit test..")

		expected_json = self.get_example_result_json("get_smiles_from_name")  # gets expected response json for mock requests.post
		expected_response = {'smiles': self.test_smiles, 'format': "smiles"}  # expected response from test function

		# Testing function with a mock of web_call function:
		with patch('qed.cts_app.cts_calcs.calculator.Calculator.web_call') as service_mock:

			# Sets web_call mock to return expected json when test function calls it:
			service_mock.return_value = expected_json

			# Calls test function, which will use the mock web_call for unit testing:
			response = self.calc_obj.get_smiles_from_name(self.test_chemical)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_response)
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_check_response_for_errors(self):
		"""
		Testing calculator module's check_response_for_errors function, which calls
		calculator function that checks for error keys in JchemWS responses.
		"""

		print(">>> Running calculator check_response_for_errors unit test..")

		test_dict = {'data': None}  # testing w/ arbitrary dict
		expected_response = {'valid': True, 'error': None}  # expected result: no error keys found

		response = self.calc_obj.check_response_for_errors(test_dict)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_response)
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_handle_error_message(self):
		"""
		Testing calculator handle_error_message function, which
		checks error keys and returns an error message for the user.
		"""

		print(">>> Running calculator handle_error_message unit test..")

		test_dict = {'errorCode': 3}  # example error key from jchemws
		expected_response = "Chemical cannot be standardized.."  # expected error message given test error code

		response = self.calc_obj.handle_error_messages(test_dict)

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_response)
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_web_call(self):
		"""
		Testing calculator module's web_call function, which calls
		JchemWS for all requests.
		"""

		print(">>> Running calculator web_call unit test..")

		mock_object = {'test': True}  # mock response from requests.post in test function
		expected_response = {'test': True}  # expected response from test function

		# Testing function with a mock of requests.post, which is used by test function:
		with patch('qed.cts_app.cts_calcs.calculator.requests.post') as service_mock:

			# Sets requests.post mock to return expected response content when test function calls it:
			service_mock.return_value.content = json.dumps(mock_object)

			# Calls test function, which will use the mock requests.post for unit testing:
			response = self.calc_obj.web_call("/fake/url", mock_object)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_response)
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_nodeWrapper(self):
		"""
		Testing calculator module's nodeWrapper function, which calls
		JchemWS for a chemical's mass.
		"""

		print(">>> Running calculator nodeWrapper unit test..")

		# Expected HTML to be retured from nodeWrapper:
		expected_response = '<div style="background-color:white;">iVBORw0KGgoAAAANSUhEUgAAAZAAAADzCAYAAACoqECMAAAQK0lEQVR42u3db2jV52LA8UON/2hweSGdA8uV4YtSZMSRoRRluRBBbkVCptSxFrJOii984QuHssgIyOpAuBWk+MIXMvpCaEt9kdFeKCwr4rwl5fZSuznmIC/sndCM6zR3jTRuvz1P8py7p8fk5HfMH3PO+Qgfojk5yTlPjr9vfr/n9/xSKYqiAgCNMggACAgAAgKAgAAgIAAgIAAICAACAoCAACAgACAgAAgIAAICgIAAICAAICAACAgAAgKAgAAgIAAgIAAICAACAoCAACAgACAgAAgIAAICgIAAICDQsL949dVXb6xfv/7Rxo0bv3vzzTd/Ft73Rz94gc28xOZ58dW5DRAQWtS77777dykY3cFzQUfw5y+++OKvwttdAgICAnMZiHse89z2Rm9v7xcCAgICT3jrrbc+Dm8PznN73BvZLiAgIPCETZs2PQxv15V6gQkICAhUrV279vvSL7DwEqvHeIKAYA/EHggICNQ3NDT0YXh7oM7HHGgwIL8XP2c8FTgaHh5+P7zvd401CAitp39gYOCzeW7bfuHChfcaCci+ffs+D2+PpVOB4yT80OHDh0eNMwgILejUqVMfReHvL6f3bYoR6Onp+Tq8fXGxh7AamWcBBITmc3Lnzp3/smbNmsfx0FNcXBje96NFzoHEPZA/TYexjDEICCzs7bffvvr888//JsYoX80OCAiUMj4+/re7du36yliAgECjOsyBgIDA01j3wgsvTBgHEBCoK07Eh7f7q5Po33777d9ExgYEBBby+2+88can8feKxEn0OAdiTEBAABAQABAQlvJFVKlcDEaMBQgINBqQ0WDcWICAgIAAAoKAAAKCgAACgoAAAoIXkYCAgICAAAKCgAACgoAAAoKAAAICAgICAgICCAgCAggIAgIICAICCAgICCAgCAggIAgIICAICCAgtENAbgTfGgsQECgbjq7gdPA4+N+g07iAgEC9cHQEJ4P7QXzHfwX3jA0ICNSLR39wJ4Xjevr3P5oDAQGB+cJxILidwjEeHIl7Ium2m8F/GCcQEMjD0ZP2NOI/JoMztXMd4c9/ptuvBd3GDQSE9g7HjuCTYDqYCs7GSfN5PnZP8LP0sfEdfx9sM44gILRXODYHl7IYXA22lrzvoWxv5HH6PF3GFQSE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</div>'

		# Test inputs for nodeWrapper's smilesToImage function call:
		test_input = self.node_wrapper_test_input

		# Mock response for smilesToImage function called within test function
		mock_json_smiles_to_image = self.get_example_result_json('smiles_to_image')

		with patch('qed.cts_app.cts_calcs.calculator.Calculator.smilesToImage') as service_mock:

			# Sets smilesToImage mock to return expected json when test function calls it:
			service_mock.return_value = mock_json_smiles_to_image

			# Calls test function, which will use the mock smilesToImage for unit testing:
			response = self.calc_obj.nodeWrapper(self.test_chemical, test_input['height'],
				test_input['width'], test_input['scale'], test_input['key'], test_input['img_type'])

		try:
			# Compares function response with expected json:
			self.assertEqual(response, expected_response)
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))



	def test_popupBuilder(self):
		"""
		Testing calculator module's popupBuilder function, which calls
		JchemWS for a chemical's mass.
		"""

		print(">>> Running calculator popupBuilder unit test..")

		# Expected HTML to be retured from mock nodeWrapper function, which is called by test function:
		mock_html = '<div style="background-color:white;">iVBORw0KGgoAAAANSUhEUgAAAZAAAADzCAYAAACoqECMAAAQK0lEQVR42u3db2jV52LA8UON/2hweSGdA8uV4YtSZMSRoRRluRBBbkVCptSxFrJOii984QuHssgIyOpAuBWk+MIXMvpCaEt9kdFeKCwr4rwl5fZSuznmIC/sndCM6zR3jTRuvz1P8py7p8fk5HfMH3PO+Qgfojk5yTlPjr9vfr/n9/xSKYqiAgCNMggACAgAAgKAgAAgIAAgIAAICAACAoCAACAgACAgAAgIAAICgIAAICAAICAACAgAAgKAgAAgIAAgIAAICAACAoCAACAgACAgAAgIAAICgIAAICDQsL949dVXb6xfv/7Rxo0bv3vzzTd/Ft73Rz94gc28xOZ58dW5DRAQWtS77777dykY3cFzQUfw5y+++OKvwttdAgICAnMZiHse89z2Rm9v7xcCAgICT3jrrbc+Dm8PznN73BvZLiAgIPCETZs2PQxv15V6gQkICAhUrV279vvSL7DwEqvHeIKAYA/EHggICNQ3NDT0YXh7oM7HHGgwIL8XP2c8FTgaHh5+P7zvd401CAitp39gYOCzeW7bfuHChfcaCci+ffs+D2+PpVOB4yT80OHDh0eNMwgILejUqVMfReHvL6f3bYoR6Onp+Tq8fXGxh7AamWcBBITmc3Lnzp3/smbNmsfx0FNcXBje96NFzoHEPZA/TYexjDEICCzs7bffvvr888//JsYoX80OCAiUMj4+/re7du36yliAgECjOsyBgIDA01j3wgsvTBgHEBCoK07Eh7f7q5Po33777d9ExgYEBBby+2+88can8feKxEn0OAdiTEBAABAQABAQlvJFVKlcDEaMBQgINBqQ0WDcWICAgIAAAoKAAAKCgAACgoAAAoIXkYCAgICAAAKCgAACgoAAAoKAAAICAgICAgICCAgCAggIAgIICAICCAgICCAgCAggIAgIICAICCAgtENAbgTfGgsQECgbjq7gdPA4+N+g07iAgEC9cHQEJ4P7QXzHfwX3jA0ICNSLR39wJ4Xjevr3P5oDAQGB+cJxILidwjEeHIl7Ium2m8F/GCcQEMjD0ZP2NOI/JoMztXMd4c9/ptuvBd3GDQSE9g7HjuCTYDqYCs7GSfN5PnZP8LP0sfEdfx9sM44gILRXODYHl7IYXA22lrzvoWxv5HH6PF3GFQSE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</div>'

		# Test input for popupBuilder function call:
		test_dict = self.popup_builder_test_input

		# Expected response from test function:
		expected_response = self.popup_builder_test_input
		expected_html = self.get_example_result_html("popup_builder")  # expected html
		expected_response['html'] = expected_html

		with patch('qed.cts_app.cts_calcs.calculator.Calculator.nodeWrapper') as service_mock:

			# Sets nodeWrapper mock to return expected json when test function calls it:
			service_mock.return_value = mock_html

			# Calls test function, which will use the mock nodeWrapper for unit testing:
			response = self.calc_obj.popupBuilder(test_dict, list(test_dict.keys()))

		try:
			# Compares function response with expected json:
			# TODO: Assert full objects are equal, currently having trouble getting the expected HTML right..
			self.assertEqual(response['smiles'], expected_response['smiles'])
		finally:
			tab = [[response], [expected_response]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))