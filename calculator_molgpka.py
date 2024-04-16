import requests
import json
import logging
import os

from .calculator import Calculator


class MolgpkaCalc(Calculator):
    """
    Handles requests to and from cts-pkasolver.
    """

    def __init__(self):
        Calculator.__init__(self)
        self.baseUrl = os.environ.get('CTS_MOLGPKA_SERVER', 'http://localhost')
        self.urlStruct = "/molgpka/data"
        self.molgpka_api_url = self.baseUrl + self.urlStruct
        self.name = 'molgpka'
        self.timeout = 10  # request timeout in seconds
        self.meta_info = {
            "metaInfo": {
                "model": "molgpka",
                "collection": "qed",
                "modelVersion": "1.1",
                "description": "Fast and accurate prediction of the pKa values of small molecules is important in the drug discovery process since the ionization state of a drug has significant influence on its activity and ADME-Tox properties. MolGpKa is a tool for pKa prediction using graph-convolutional neural network model. The model works by learning pKa related chemical patterns automatically and building reliable predictors with learned features.",
                "timestamp": self.gen_jid(),
                "url": "https://github.com/Xundrug/MolGpKa"
            }
        }
        self.response_obj = {
            'calc': "molgpka",  # todo: change to metabolizer, change in template too
            'prop': "pchem",
            'data': None,
            'chemical': None,
            'request_post': None            
        }

    def validate_response(self, results):
        """
        Validate response content from molgpka.
        """
        if results.get("status") != True:
            # TODO: Send error like other calculator_*.py modules.
            error_response = {
                "error": "Error making request to cts-pkasolver",
                "valid": False
            }
            return error_response

        return results

    def data_request_handler(self, request_dict):

        chemical = request_dict["chemical"]

        # TODO: Any sort of SMILES validation??s

        _response_obj = dict(self.response_obj)

        post_data = {
            "smiles": chemical
        }
        
        results = None

        try:
            response = requests.get(self.molgpka_api_url, params=post_data, timeout=self.timeout)
            results = json.loads(response.content)
            results = self.validate_response(results)
        except Exception as e:
            logging.warning("calculator_molgpka exception: {}".format(e))
            _response_obj.update({"valid": False, "error": "Error getting data from molgpka"})
            return _response_obj

        _response_obj['data'] = results
        _response_obj['chemical'] = request_dict.get('chemical')
        _response_obj['request_post'] = request_dict

        return _response_obj
