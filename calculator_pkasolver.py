import requests
import json
import logging
import os

from .calculator import Calculator


class PkaSolverCalc(Calculator):
    """
    Handles requests to and from cts-pkasolver.
    """

    def __init__(self):
        Calculator.__init__(self)
        self.baseUrl = os.environ.get('CTS_PKASOLVER_SERVER', "http://localhost:8181")
        self.urlStruct = "/pkasolver/data"
        self.pkasolver_api_url = self.baseUrl + self.urlStruct
        self.name = 'pkasolver'
        self.timeout = 600  # request timeout in seconds
        self.products_list = []  # single-level list of products with genKey
        self.metID = 0  # unique id for each node
        self.meta_info = {
            "metaInfo": {
                "model": "pkasolver",
                "collection": "qed",
                "modelVersion": "0.3",
                "description": "pKasolver is a package that enables prediction of pKa values of small molecules via graph convolutional networks.",
                "timestamp": self.gen_jid(),
                "url": "https://github.com/mayrf/pkasolver"
            }
        }
        self.response_obj = {
            'calc': "pkasolver",  # todo: change to metabolizer, change in template too
            'prop': "pchem",
            # 'node': None,
            'data': None,
            # 'total_products': None,
            'chemical': None,
            # 'workflow': 'gentrans',
            # 'run_type': 'batch',
            'request_post': None            
        }

    def validate_response(self, response):
        """
        Validate response content from cts_pkasolver.
        Example repsonse: {
            "pka_list": [9.63],
            "status": true
        }
        """
        if response.get("status") != True:
            # TODO: Send error like other calculator_*.py modules.
            error_response = {
                "data": "Error making request to cts-pkasolver",
                "valid": False
            }
            return error_response

        # Gets pka values as a list:
        pka_list = response.get("pka_list")

        if not pka_list or len(pka_list) < 1:
            error_response = {
                "data": "Error making request to cts-pkasolver",
                "valid": False
            }
            return error_response

        return pka_list

    def data_request_handler(self, request_dict):

        chemical = request_dict["chemical"]
        # data_type = request_dict.get("data_type")

        # TODO: Any sort of SMILES validation??

        _response_obj = dict(self.response_obj)

        post_data = {
            "smiles": chemical
        }


        
        # TODO: Error handling for request. If pkasolver fails continue on with the
        # normal speciation output.



        try:
            response = requests.get(self.pkasolver_api_url, params=post_data, timeout=self.timeout)
        except requests.exceptions.RequestException as e:
            logging.warning("calculator_pkasolver exception: {}".format(e))
            _response_obj.update({'error': "Error getting data from pkasolver"})
            return _response_obj

        results = json.loads(response.content)

        microspecies = []

        # Get chem info for returned microspecies:
        for key, smiles in results.get("species", {}).items():
            chem_info = self.getStructInfo(smiles)
            ms_key = "microspecies" + str(int(key) + 1)
            chem_info.update({"key": ms_key})
            microspecies.append(chem_info)

        results["pka_microspecies"] = microspecies
        del results["species"]

        _response_obj['data'] = results
        _response_obj['chemical'] = request_dict.get('chemical')
        _response_obj['request_post'] = request_dict

        return _response_obj
