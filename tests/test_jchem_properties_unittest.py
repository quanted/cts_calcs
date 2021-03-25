"""
Unit testing for jchem_properties module, which
handles JchemWS property requests.
"""

import unittest
import json
import os
import inspect
import datetime
import logging
import sys
from tabulate import tabulate

_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(
    1, os.path.join(_path, "..", "..", "..", "..")
)  # adds qed project to sys.path

# local requirements (running pytest at qed level):
if 'cts_celery' in _path:
	from qed.cts_celery.cts_calcs.jchem_properties import JchemProperty
elif 'cts_app' in _path:
	from qed.cts_app.cts_calcs.jchem_properties import JchemProperty

from qed.temp_config.set_environment import DeployEnv



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
		self.filename_structure = "mock_json/jchem_properties_result_{}.json"  # {} = property (e.g., pka, tautomer)

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
			assert math.isclose(acidic_vals, acidic_expected_results, rel_tol=1e-6, abs_tol=0)
			assert math.isclose(basic_vals, basic_expected_results, rel_tol=1e-6, abs_tol=0)
		finally:
			tab = [result, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_pka_parent(self):
		"""
		Testing pKa class getParent function, which grabs
		the parent chemical's image from JchemWS in base64 encoding.
		"""
		print(">>> Running pka parent test..")

		expected_smiles = ["CC(N)C(O)=O"]

		test_json = self.get_example_result_json('pka')
		test_obj = self.jc.getPropObject('pKa')
		test_obj.results = test_json

		try:
			parent_image = test_obj.getParent(True)
			print("pka function getParent() test passed!")
		finally:
			print("\n")
			
		return



	def test_pka_microspecies(self):
		"""
		Testing pKa class getMicrospecies function, which
		gets microspecies list of images and chem info.
		"""
		print(">>> Running pka microspecies test..")

		test_json = self.get_example_result_json('pka')
		test_obj = self.jc.getPropObject('pKa')
		test_obj.results = test_json

		try:
			ms_image = test_obj.getMicrospecies(True)
			print("pka function getMicrospecies() test passed!")
		finally:
			print("\n")
			
		return



	def test_pka_chart(self):
		"""
		Testing pKa class getChartData, which builds
		a dictionary of microspecies data across pH for
		the speciation workflow.
		"""
		print(">>> Running pka chart test..")

		expected_ms_total = 3  # expected number of species from pka test result json
		expected_ms1_num = 141  # expected number of values for microspecies1 plot
		expected_ms2_num = 141  # expected number of values for microspecies2 plot
		expected_ms3_num = 141  # expected number of values for microspecies3 plot
		expected_results = [expected_ms_total, expected_ms1_num, expected_ms2_num, expected_ms3_num]

		test_json = self.get_example_result_json('pka')
		test_obj = self.jc.getPropObject('pKa')
		test_obj.results = test_json

		chart_data = test_obj.getChartData()

		ms_total = len(chart_data)
		ms1_num = len(chart_data['microspecies1'])
		ms2_num = len(chart_data['microspecies2'])
		ms3_num = len(chart_data['microspecies3'])
		results = [ms_total, ms1_num, ms2_num, ms3_num]

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-2, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_isoelectricpoint_get_isopt(self):
		"""
		Testing IsoelectricPoint class getIsoelectricPoint function,
		which grabs the float value.
		"""
		print(">>> Running isoelectricpoint getip test..")

		expected_results = [5.975982685926768]

		test_json = self.get_example_result_json('isopt')
		test_obj = self.jc.getPropObject('isoelectricPoint')
		test_obj.results = test_json

		isopt = test_obj.getIsoelectricPoint()

		results = [isopt]

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-12, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_isoelectricpoint_chart(self):
		"""
		Testing IsoelectricPoint class getChartData function,
		which grabs isoelectricpoint chart data
		"""
		print(">>> Running isoelectricpoint getip test..")

		expected_results = [141]  # number of x,y pairs for isoelectric point plot

		test_json = self.get_example_result_json('isopt')
		test_obj = self.jc.getPropObject('isoelectricPoint')
		test_obj.results = test_json

		chart_data = test_obj.getChartData()

		results = [len(chart_data)]  # number of x,y pairs for isoelectric point

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-2, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_major_microspecies(self):
		"""
		Testing major microspecies class getMajorMicrospecies function,
		which grabs images of the major microspecies
		"""
		print(">>> Running major microspecies test..")

		# major microspecies for test chemical, and base64 encoded image of major ms:
		expected_results = ["CC([NH3+])C([O-])=O"]

		test_json = self.get_example_result_json('majormicrospecies')
		test_obj = self.jc.getPropObject('majorMicrospecies')
		test_obj.results = test_json

		major_ms_smiles = test_obj.results['result']['structureData']['structure']
		major_ms_image = test_obj.getMajorMicrospecies(True)

		results = [major_ms_smiles]

		try:
			for i in range(0, len(expected_results)):
				self.assertEqual(expected_results[i], results[i])
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_tautomerization(self):
		"""
		Testing Tautomerization class getTautomers function,
		which grabs tautomers of parent chemical.
		"""
		print(">>> Running tautomers test..")

		expected_results = [2, 99.66, 0.33]

		test_json = self.get_example_result_json('tautomerization')
		test_obj = self.jc.getPropObject('tautomerization')
		test_obj.results = test_json

		tauts = test_obj.getTautomers(True)

		taut_dist_1 = tauts[0]['dist']
		taut_dist_2 = tauts[1]['dist']

		results = [len(tauts), taut_dist_1, taut_dist_2]

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-2, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_stereoisomer(self):
		"""
		Testing stereoisomer class getStereoisomers function,
		which grabs stereoisomers of parent chemical.
		"""
		print(">>> Running stereoisomers test..")

		expected_results = [2, "C[C@H](N)C(O)=O", "C[C@@H](N)C(O)=O"]

		test_json = self.get_example_result_json('stereoisomer')
		test_obj = self.jc.getPropObject('stereoisomer')
		test_obj.results = test_json

		stereos = test_obj.getStereoisomers(True)

		stereo_1 = test_obj.results['result'][0]['structureData']['structure']
		stereo_2 = test_obj.results['result'][1]['structureData']['structure']

		results = [len(stereos), stereo_1, stereo_2]

		try:
			for i in range(0, len(expected_results)):
				self.assertEqual(expected_results[i], results[i])
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_solubility(self):
		"""
		Testing solubility class functions.
		"""
		print(">>> Running solubility tests..")

		expected_results = [1000*169.23086532148588, 0.2801]  # intrinsic sol and solubility at pH 7.0

		test_json = self.get_example_result_json('solubility')
		test_obj = self.jc.getPropObject('solubility')
		test_obj.results = test_json

		intrinsic_sol = test_obj.getIntrinsicSolubility()  # gets intrinsic solubility
		ph_sol = test_obj.getPHDependentSolubility(7.0)  # gets solubility at pH

		results = [intrinsic_sol, ph_sol]

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-6, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_logp(self):
		"""
		Testing LogP class functions.
		"""
		print(">>> Running logp tests..")

		expected_results = [-0.5787113146666667]

		test_json = self.get_example_result_json('logp')
		test_obj = self.jc.getPropObject('logP')
		test_obj.results = test_json

		logp = test_obj.getLogP()

		results = [logp]

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-12, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



	def test_logd(self):
		"""
		Testing LogD class functions.
		"""
		print(">>> Running logd tests..")

		expected_results = []

		test_json = self.get_example_result_json('logd')
		test_obj = self.jc.getPropObject('logD')
		test_obj.results = test_json

		# logp = test_obj.getLogP()
		logd = test_obj.getLogD(7.0)

		results = []

		try:
			assert math.isclose(results, expected_results, rel_tol=1e-12, abs_tol=0)
		finally:
			tab = [results, expected_results]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))
			
		return



if __name__ == '__main__':
	unittest.main()