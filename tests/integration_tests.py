import sys
import os
import requests
import unittest
import logging
import json
# import numpy.testing as npt
# import linkcheck_helper
# from . import linkcheck_helper
# from .test_objects import get_post_object

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(1, os.path.join(PROJECT_ROOT, "..", ".."))

from cts_calcs.tests.test_objects import get_post_object
# from django.test import Client, TestCase

from temp_config.set_environment import DeployEnv

# Determine env vars to use:
runtime_env = DeployEnv()
runtime_env.load_deployment_environment()

# import django
# django.setup()

# servers = [os.getenv("CTS_REST_SERVER")]
servers = ["https://qedlinux1stg.aws.epa.gov"]

print("SERVERS: {}".format(servers))

pages = [
    "/cts/",
    "/cts/rest/",
    "/cts/chemspec",
    "/cts/pchemprop",
    "/cts/gentrans",
    "/cts/chemspec/input",
    "/cts/pchemprop/input",
    "/cts/gentrans/input",
    "/cts/chemspec/batch",
    "/cts/pchemprop/batch",
    "/cts/gentrans/batch"
]

workflow_endpoints = [
    "/cts/chemspec/output/"
]

api_endpoints = [
    "/cts/rest/chemaxon/run",
    "/cts/rest/epi/run",
    "/cts/rest/measured/run",
    "/cts/rest/testws/run",
    "/cts/rest/opera/run"
]

api_endpoints_map = {
    "chemaxon": "/cts/rest/chemaxon/run",
    "epi": "/cts/rest/epi/run",
    "measured": "/cts/rest/measured/run",
    "testws": "/cts/rest/testws/run",
    "opera": "/cts/rest/opera/run",
    "metabolizer": "/cts/rest/metabolizer/run"
}

# following are lists of url's to be processed with tests below
check_pages = [s + p for s in servers for p in pages]

print("checking pages: {}".format(check_pages))
# print("checking workflow outputs at: {}".format(workflow_test_urls))


class TestCTSPages(unittest.TestCase):
# class TestCTSPages(TestCase):
    """
    this testing routine accepts a list of pages and performs a series of unit tests that ensure
    that the web pages are up and operational on the server.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_cts_200(self):
        """
        Tests basic HTML pages.
        """
        test_name = "Checks page access "
        responses = [requests.get(p).status_code for p in check_pages]
        for response in responses:
            self.assertEqual(response, 200, "200 error")
        # npt.assert_array_equal(response, 200, '200 error', True)

    # def test_cts_workflows():
    #     """
    #     Tests workflow outputs.
    #     """

        # test_client = Client()
    #     test_name = "Checks workflow outputs "
    #     response = []

    #     for p in workflow_endpoints:
    #         workflow = p.split('/')[2]  # assuming /cts/model/output url structure, grabbing "model" part
    #         url = os.getenv('CTS_REST_SERVER') + p
    #         post_data = get_post_object(workflow)

    # print("~~~ {}".format(test_name))
    #         print("url: {}".format(url))
    #         print("post: {}".format(post_data))

    #         res = requests.post(url, data=post_data)

    #         self.assertEqual(res.status_code, 200, "200 error")
            # response.append(res.status_code)
        # npt.assert_array_equal(response, 200, '200 error', True)

    def test_cts_metabolizer_endpoint(self):
        """
        Tests Metabolizer API endpoint.
        """
        test_name = "Checks CTS API Metabolizer endpoint "
        # test_client = Client()
        response = []
        calc_name = "metabolizer"

        calc_endpoint = api_endpoints_map[calc_name]
        url = servers[0] + calc_endpoint
        post_data = get_post_object(calc_name)

        print("~~~ {}".format(test_name))
        print("url: {}".format(url))
        print("post: {}".format(post_data))

        res = requests.post(url, data=post_data)
        result = json.loads(res.content)

    def test_cts_api_chemaxon_endpoint(self):
        """
        Tests ChemAxon API endpoint.
        """
        test_name = "Checks CTS API ChemAxon endpoint "
        # test_client = Client()
        response = []
        calc_name = "chemaxon"

        calc_endpoint = api_endpoints_map[calc_name]
        url = servers[0] + calc_endpoint
        post_data = get_post_object(calc_name)

        print("~~~ {}".format(test_name))
        print("url: {}".format(url))
        print("post: {}".format(post_data))

        res = requests.post(url, data=post_data)
        result = json.loads(res.content)

        # npt.assert_equal(False, 'error' in result, verbose=True)
        self.assertEqual(False, "error" in result)

    def test_cts_api_epi_endpoint(self):
        """
        Tests EPI API endpoint.
        """
        test_name = "Checks CTS API EPI Suite endpoint "
        # test_client = Client()
        response = []
        calc_name = "epi"

        calc_endpoint = api_endpoints_map[calc_name]
        url = servers[0] + calc_endpoint
        post_data = get_post_object(calc_name)

        print("~~~ {}".format(test_name))
        print("url: {}".format(url))
        print("post: {}".format(post_data))

        res = requests.post(url, data=post_data)
        result = json.loads(res.content)

        # npt.assert_equal(False, 'error' in result, verbose=True)
        self.assertEqual(False, 'error' in result)

    def test_cts_api_testws_endpoint(self):
        """
        Tests TESTWS API endpoint.
        """
        test_name = "Checks CTS API TESTWS endpoint "
        # test_client = Client()
        response = []
        calc_name = "testws"

        calc_endpoint = api_endpoints_map[calc_name]
        url = servers[0] + calc_endpoint
        post_data = get_post_object(calc_name)

        print("~~~ {}".format(test_name))
        print("url: {}".format(url))
        print("post: {}".format(post_data))

        res = requests.post(url, data=post_data)
        result = json.loads(res.content)

        # npt.assert_equal(False, 'error' in result, verbose=True)
        self.assertEqual(False, 'error' in result)

    def test_cts_api_opera_endpoint(self):
        """
        Tests OPERA API endpoint.
        """
        test_name = "Checks CTS API OPERA endpoint "
        # test_client = Client()
        response = []
        calc_name = "opera"

        calc_endpoint = api_endpoints_map[calc_name]
        url = servers[0] + calc_endpoint
        post_data = get_post_object(calc_name)

        print("~~~ {}".format(test_name))
        print("url: {}".format(url))
        print("post: {}".format(post_data))

        res = requests.post(url, data=post_data)
        result = json.loads(res.content)

        # npt.assert_equal(True, result.get('valid'), verbose=True)
        self.assertEqual(True, result.get('valid'))


    def test_cts_api_measured_endpoint(self):
        """
        Tests Measured API endpoint.
        """
        test_name = "Checks CTS API EPI's Measured endpoint "
        # test_client = Client()
        response = []
        calc_name = "measured"

        calc_endpoint = api_endpoints_map[calc_name]
        url = servers[0] + calc_endpoint
        post_data = get_post_object(calc_name)

        print("~~~ {}".format(test_name))
        print("url: {}".format(url))
        print("post: {}".format(post_data))

        res = requests.post(url, data=post_data)
        result = json.loads(res.content)

        # npt.assert_equal(False, 'error' in result, verbose=True)
        self.assertEqual(False, 'error' in result)



# unittest will
# 1) call the setup method,
# 2) then call every method starting with "test",
# 3) then the teardown method
if __name__ == '__main__':
    unittest.main()