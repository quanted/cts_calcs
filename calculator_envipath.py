import requests
import json
import logging
import os

from .calculator import Calculator


class EnvipathCalc(Calculator):
    """
    Handles requests to and from cts-envipath.
    """

    def __init__(self):
        Calculator.__init__(self)
        self.baseUrl = os.environ.get('CTS_ENVIPATH_SERVER', "http://localhost:5003")
        self.urlStruct = "/envipath/rest/run"
        self.envipath_api_url = self.baseUrl + self.urlStruct
        self.name = 'envipath'
        self.timeout = 600  # request timeout in seconds
        self.products_list = []  # single-level list of products with genKey
        self.metID = 0  # unique id for each node
        self.meta_info = {
            "metaInfo": {
                "model": "envipath",
                "collection": "qed",
                "modelVersion": "",
                "description": "enviPath is a database and prediction system for the microbial biotransformation of organic environmental contaminants. The database provides the possibility to store and view experimentally observed biotransformation pathways. The pathway prediction system provides different relative reasoning models to predict likely biotransformation pathways and products.",
                "timestamp": self.gen_jid(),
                "url": "https://envipath.org"
            }
        }
        self.response_obj = {
            'calc': "envipath",  # todo: change to metabolizer, change in template too
            'prop': "products",
            'node': None,
            'data': None,
            'total_products': None,
            'chemical': None,
            'workflow': 'gentrans',
            'run_type': 'batch',
            'request_post': None            
        }

    def recursive(self, jsonDict, gen_limit, unranked=False):
        """
        Starting point for walking through
        metabolites dictionary and building json
        that thejit (visualization javascript
        library) understands
        """
        root = jsonDict['data']
        reDict = {}
        reDict.update({
            'tree': self.traverse(root, gen_limit, unranked),
            'total_products': self.metID - 1  # subtract out the parent for "total products" value
        })
        self.metID = 0  # resets the metID attribute
        return json.dumps(reDict)


    def traverse(self, root, gen_limit, unranked=False):
        """
        For gentrans model output - products tree
        Uses JIT spacetree library on output page; this
        function walks the metabolizer response tree and uses
        keys for plotting a JIT spacetree.
        (https://philogb.github.io/jit/)

        id - generation ID (e.g., 2.4.3 - parent 2's 4th child's 3rd child..)
        name - product image (node) on the spacetree
        data - shows up when the house hovers over a product
        """

        self.metID += 1
        _products_dict = {}

        logging.info("metabolites: {}".format(self.metID))

        if self.metID == 1:

            _parent = root['smiles']

            _products_dict.update({
                "id": self.metID,
                "name": "<img class='blank_node' src='" + settings.STATIC_URL + "cts_app/images/loader_node.gif' />",
                "data": {
                    'smiles': _parent,
                    'routes': "",
                    'accumulation': "N/A" if unranked else round(root.get('accumulation'), 4),
                    'production': "N/A" if unranked else round(root.get('production'), 4),
                    'globalAccumulation': "N/A" if unranked else round(root.get('globalAccumulation'), 4),
                    'likelihood': "N/A" if unranked else root.get('likelihood')
                },
                "children": []
            })
            
        else:
            _products_dict.update({
                "id": self.metID,
                "name": "<img class='blank_node' src='" + settings.STATIC_URL + "cts_app/images/loader_node.gif' />",
                "data": {
                    'smiles': root['smiles'],
                    'routes': "",
                    'accumulation': "N/A" if unranked else round(root.get('accumulation'), 4),
                    'production': "N/A" if unranked else round(root.get('production'), 4),
                    'globalAccumulation': "N/A" if unranked else round(root.get('globalAccumulation'), 4),
                    'likelihood': "N/A" if unranked else root.get('likelihood')
                },
                "children": []
            })

        if 'metabolites' in root and len(root['metabolites']) > 0:
            for metabolite in root['metabolites']:
                root2 = metabolite
                _products_dict['children'].append(self.traverse(root2, gen_limit, unranked))

        return _products_dict



    def data_request_handler(self, request_dict):

        # metabolizer_data = request_dict.get("metabolizer_post")
        chemical = request_dict["chemical"]
        # prop = metabolizer_data.get("prop")
        gen_limit = int(request_dict.get("gen_limit", 1))

        _response_obj = dict(self.response_obj)

        post_data = {
            "smiles": chemical,
            "gen_limit": gen_limit
        }

        try:
            response = requests.post(self.envipath_api_url, json=post_data, timeout=self.timeout)
        except requests.exceptions.RequestException as e:
            logging.warning("calculator_envipath exception: {}".format(e))
            _response_obj.update({'error': "Error getting data from Envipath"})
            return _response_obj

        _results = self.recursive(json.loads(response.content), int(request_dict['gen_limit']), True)

        _products_data = json.loads(_results)

        _response_obj['node'] = request_dict.get('node')
        _response_obj['data'] = _products_data['tree']
        _response_obj['total_products'] = _products_data['total_products']
        _response_obj['chemical'] = request_dict.get('chemical')
        _response_obj['request_post'] = request_dict

        return _response_obj
