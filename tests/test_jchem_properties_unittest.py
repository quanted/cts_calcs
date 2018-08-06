import unittest
import json
import os
import inspect
import datetime
import logging
import sys
from numpy import testing as npt
from tabulate import tabulate

# local requirements (running pytest at qed level):
try:
	from qed.cts_app.cts_calcs.jchem_properties import JchemProperty
	from qed.temp_config.set_environment import DeployEnv
except ModuleNotFoundError as e:
	logging.warning(">>> Exception in test_calculator_unittest.py importing local reqs: {}".format(e))
	logging.critical(">>> Ignoring test_calculator_unittest.py for now (i.e., quitting this test script)!!")
	quit()


# class TestJchemProperties(object):
class TestJchemProperties(unittest.TestCase):

	print("jchem properties unittests conducted at " + str(datetime.datetime.today()))

	def setUp(self):
		"""
		Setup routine for Kabam unit tests.
		:return:
		"""

		# Sets up runtime environment:
		runtime_env = DeployEnv()
		runtime_env.load_deployment_environment()

		# Defines filename structure for example JSON results from jchem WS:
		self.filename_structure = "jchem_unittest_object_{}.json"  # {} = property (e.g., pka, tautomer)

		self.jc = JchemProperty()  # shared instance of jchem_properties module



	def tearDown(self):
		"""
		Teardown routine for Kabam unit tests.
		:return:
		"""
		pass
		# teardown called after each test
		# e.g. maybe write test results to some text file

	

	def create_jchem_object(self, prop):
		# # create empty pandas dataframes to create empty object for testing
		# df_empty = pd.DataFrame()
		# # create an empty kabam object
		# kabam_empty = Kabam(df_empty, df_empty)
		# return kabam_empty
		pass



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



	def test_pka_acidic_basic(self):
		"""
		Testing pKa class acidic/basic pKa values. Loads .json file of sample
		Jchem WS pKa result, and grab basic and acidic values using jchem_properties pKa class.
		"""
		print(">>> Running pka acidic/basic test..")

		acidic_expected_results = [2.474897471379417]
		basic_expected_results = [9.476905287351835]

		test_json = self.get_example_result_json('pka')
		test_obj = self.jc.getPropObject('pKa')
		test_obj.results = test_json

		acidic_vals = test_obj.getMostAcidicPka()
		basic_vals = test_obj.getMostBasicPka()

		result = [acidic_vals, basic_vals]
		expected_results = [acidic_expected_results, basic_expected_results]

		try:
			npt.assert_allclose(acidic_vals, acidic_expected_results, rtol=1e-6, atol=0, err_msg='', verbose=True)
			npt.assert_allclose(basic_vals, basic_expected_results, rtol=1e-6, atol=0, err_msg='', verbose=True)
		finally:
			tab = [result, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	# NOTE: getParent() function makes request to get chem info. Write tests
	# for Calculator class before jchem_properties (jp uses calc).

	# def test_pka_parent(self):
	# 	"""
	# 	Testing pKa class, getting parent molecule. Loads .json file of sample
	# 	Jchem WS pKa result, and grabs parent chemical using jchem_properties pKa class.
	# 	"""
	# 	print(">>> Running pka parent test..")

	# 	expected_result = ["CC(N)C(O)=O"]

	# 	test_json = self.get_example_result_json('pka')
	# 	test_obj = self.jc.getPropObject('pKa')
	# 	test_obj.results = test_json

	# 	result = []

	# 	try:

	# 		npt.assert_allclose(acidic_vals, acidic_expected_results, rtol=1e-6, atol=0, err_msg='', verbose=True)
	# 		npt.assert_allclose(basic_vals, basic_expected_results, rtol=1e-6, atol=0, err_msg='', verbose=True)
	# 	finally:
	# 		tab = [result, expected_results]
	# 		print("\n")
	# 		print(inspect.currentframe().f_code.co_name)
	# 		print(tabulate(tab, headers='keys', tablefmt='rst'))
	# 	return



if __name__ == '__main__':
	unittest.main()