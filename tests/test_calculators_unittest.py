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



	# def test_tautomerization(self):
	# 	"""
	# 	Testing Tautomerization class getTautomers function,
	# 	which grabs tautomers of parent chemical.
	# 	"""
	# 	print(">>> Running tautomers test..")

	# 	expected_results = [2, 99.66, 0.33]

	# 	test_json = self.get_example_result_json('tautomerization')
	# 	test_obj = self.jc.getPropObject('tautomerization')
	# 	test_obj.results = test_json

	# 	tauts = test_obj.getTautomers(True)

	# 	taut_dist_1 = tauts[0]['dist']
	# 	taut_dist_2 = tauts[1]['dist']

	# 	results = [len(tauts), taut_dist_1, taut_dist_2]

	# 	try:
	# 		npt.assert_allclose(results, expected_results, rtol=1e-2, atol=0, err_msg='', verbose=True)
	# 	finally:
	# 		tab = [results, expected_results]
	# 		print("\n")
	# 		print(inspect.currentframe().f_code.co_name)
	# 		print(tabulate(tab, headers='keys', tablefmt='rst'))
			
	# 	return