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
	from qed.cts_celery.cts_calcs.actorws import ACTORWS
elif 'cts_app' in _path:
	from qed.cts_app.cts_calcs.actorws import ACTORWS

from qed.temp_config.set_environment import DeployEnv



class TestACTORWS(unittest.TestCase):
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

		self.actorws_obj = ACTORWS()



	def tearDown(self):
		"""
		Teardown routine for Kabam unit tests.
		:return:
		"""
		pass
		# teardown called after each test
		# e.g. maybe write test results to some text file



	def test_make_request(self):
		"""
		Testing actorws module make_request() function.
		"""

		print(">>> Running actorws make_request unit test..")

		test_url = self.actorws_obj.chemid_url  # using chem id url for testing
		test_payload = {'identifier': self.test_chemical}

		# Expected result from actorws:
		expected_json = json.loads("""
		{
		    "DataRow": {
		        "origIdentifier": "aspirin",
		        "synGsid": 20108,
		        "synType": "Approved Name",
		        "synIdentifier": "Aspirin",
		        "inChiKey": null,
		        "smiles": null,
		        "collidingGsid": null,
		        "collidingCasrn": null,
		        "collidingPreferredName": null,
		        "trimmedWhitespace": false,
		        "trimmedLeadingZeros": false,
		        "reformattedIdentifier": false,
		        "checksum": "NA",
		        "processedAs": "LOOKUP",
		        "infoMsg": "(Retrieved by synonym)",
		        "warningMsg": null
		    }
		}
		""")

		with patch('qed.cts_app.cts_calcs.actorws.requests.get') as service_mock:

			service_mock.return_value.content = json.dumps(expected_json)  # sets expected result from actorws GET request
			service_mock.return_value.status_code = 200  # test function expects 200 status code

			response = self.actorws_obj.make_request(test_url, test_payload)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_json)
		finally:
			tab = [[response], [expected_json]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_get_chemid_results(self):
		"""
		Testing actorws module get_chemid_results() function.
		"""

		print(">>> Running actorws get_chemid_results unit test..")

		# Expected result from make_request() mock, which is called by test function:
		mock_json = json.loads("""
		{
		    "DataRow": {
		        "origIdentifier": "aspirin",
		        "synGsid": 20108,
		        "synType": "Approved Name",
		        "synIdentifier": "Aspirin",
		        "inChiKey": null,
		        "smiles": null,
		        "collidingGsid": null,
		        "collidingCasrn": null,
		        "collidingPreferredName": null,
		        "trimmedWhitespace": false,
		        "trimmedLeadingZeros": false,
		        "reformattedIdentifier": false,
		        "checksum": "NA",
		        "processedAs": "LOOKUP",
		        "infoMsg": "(Retrieved by synonym)",
		        "warningMsg": null
		    }
		}
		""")

		# Expected result from test function:
		expected_result = {
			'calc': 'actorws',
			'prop': 'chemid',
			'data': {'gsid': 20108}
		}

		# Testing function with a mock of make_request, which is called by test function:
		with patch('qed.cts_app.cts_calcs.actorws.ACTORWS.make_request') as service_mock:

			service_mock.return_value = mock_json  # sets expected result from mock function call

			response = self.actorws_obj.get_chemid_results(self.test_chemical)

		try:
			# Compares function response with expected json:
			self.assertDictEqual(response, expected_result)
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return



	def test_get_dsstox_results(self):
		"""
		Testing actorws module get_dsstox_results() function.
		"""

		print(">>> Running actorws get_dsstox_results unit test..")

		test_chemical = self.test_chemical
		test_id_type = "gsid"

		# Expected result from make_request() mock, which is called by test function:
		mock_json = json.loads("""
		{
		    "DataList": {
		        "list": [
		            {
		                "casrn": "50-78-2",
		                "gsid": 20108,
		                "dsstoxSubstanceId": "DTXSID5020108",
		                "preferredName": "Aspirin",
		                "casType": "Single Compound",
		                "qcLevel": "DSSTox_Low",
		                "qcNotes": "",
		                "qcNotes2": "single chemical compound",
		                "qcNotes3": null,
		                "source": "Public",
		                "chemnoteDeprecated": null,
		                "inchiKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
		                "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)\nAuxInfo=1/1/N:13,8,7,9,6,11,4,5,2,12,1,3,10/E:(11,12)/rA:13OCOCCCCCCOCOC/rB:s1;d2;s2;d4;s5;d6;s7;s4d8;s5;s10;d11;s11;/rC:0,-6.6665,0;.7667,-5.3332,0;0,-3.9999,0;2.3066,-5.3332,0;3.0799,-3.9999,0;4.6199,-3.9999,0;5.3866,-5.3332,0;4.6199,-6.6665,0;3.0799,-6.6665,0;2.3066,-2.6665,0;3.0799,-1.3333,0;2.3066,0,0;4.6199,-1.3333,0;\n",
		                "smiles": "CC(=O)OC1=C(C=CC=C1)C(O)=O",
		                "iupac": "C9H8O4",
		                "molFormula": "C9H8O4",
		                "molWeight": "180.159"
		            }
		        ],
		        "count": 1
		    }
		}
		""", strict=False)

		# Expected result from test function:
		expected_result = {
			'calc': 'actorws',
			'prop': 'dsstox',
			'data': {
				'gsid': 20108,
				'casrn': '50-78-2',
				'dsstoxSubstanceId': 'DTXSID5020108',
				'preferredName': 'Aspirin',
				'smiles': 'CC(=O)OC1=C(C=CC=C1)C(O)=O',
				'iupac': 'C9H8O4'
			}
		}

		# Testing function with a mock of make_request, which is called by test function:
		with patch('qed.cts_app.cts_calcs.actorws.ACTORWS.make_request') as service_mock:

			service_mock.return_value = mock_json  # sets expected result from mock function call

			response = self.actorws_obj.get_dsstox_results(test_chemical, test_id_type)

		try:
			# Compares function response with expected json:
			self.assertEqual(response['data']['dsstoxSubstanceId'], expected_result['data']['dsstoxSubstanceId'])
		finally:
			tab = [[response], [expected_result]]
			print("\n")
			print(inspect.currentframe().f_code.co_name)
			print(tabulate(tab, headers='keys', tablefmt='rst'))

		return