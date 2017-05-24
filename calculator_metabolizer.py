import requests
import json
import logging
import os
import redis

import smilesfilter
from calculator import Calculator


class MetabolizerCalc(Calculator):
    """
    Calc class for CTSWS (formerly EFS) Metabolizer.

    This was originally handled by JchemCalc and data_walks,
    but now this'll be the replacement for data_walks.

    Also put the spacetree stuff here, like image size, keys, etc.
    """

    def __init__(self, prop_name=None):
        
        Calculator.__init__(self)  # inherit Calculator base class

        # self.propsList = ['pKa', 'isoelectricPoint', 'majorMicrospecies', 'tautomerization', 'stereoisomer']  # speciation props
        
        self.baseUrl = os.environ['CTS_EFS_SERVER']
        self.name = ''
        self.max_retries = 3  # request retries if error
        self.timeout = 20  # request timeout in seconds
        self.structure = ''  # cas, smiles, iupac, etc.
        self.postData = {}
        self.results = ''
        self.methods = ['KLOP', 'VG', 'PHYS']  # kow_no_ph and kow_wph only
        self.props = ['water_sol', 'ion_con', 'kow_no_ph', 'kow_wph', 'water_sol_ph']  # available pchem props
        self.prop_name = prop_name  # prop name for MetabolizerCalc instance
        self.products_list = []  # single-level list of products with genKey

        self.metID = 0  # unique id for each node
        self.metabolite_keys = ['smiles', 'formula', 'iupac', 'mass', 'accumulation', 'production', 'transmissivity', 'generation', 'routes', 'exactMass']
        self.image_scale = 50

        self.tree_image_height = 114  # height of molecule in gentrans spacetree
        self.tree_image_width = 100  # width of molecule in gentrans spacetree

        self.gen_limit = 2

        # CTSWS Transformation Products Request
        self.transformation_request = {
            'structure': None,
            'generationLimit': 1,  # make sure to get this from front end
            'populationLimit': 0,
            'likelyLimit': 0.001,
            'transformationLibraries': ["hydrolysis", "abiotic_reduction"],  # NOTE: no transformationLibraries key:val for mammalian metabolism
            'excludeCondition': ""
        }


    def recursive(self, jsonDict, gen_limit):
        """
        Starting point for walking through
        metabolites dictionary and building json
        that thejit (visualization javascript
        library) understands
        """
        root = jsonDict['results']

        reDict = {}
        reDict.update({
            'tree': self.traverse(root, gen_limit),
            'total_products': self.metID
        })


        # Need to hit SMILES filter, then retrieve molecular info,
        # node image, and popup image for products...
        # logging.info("PRODUCTS LIST: {}".format(self.products_list))

        
        # Could make one big request to get node images from list of 
        # nodes here! But, the tree would need to be walked again so it
        # can have edit the 'name' key with the node images...


        return json.dumps(reDict)


    def traverse(self, root, gen_limit):
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
            _parent = root.keys()[0]  # start with parent metabolite

            # organizing data for spacetree
            # _products_dict.update({
            #     "id": self.metID,
            #     "name": self.nodeWrapper(_parent, self.tree_image_height, self.tree_image_width, self.image_scale, self.metID,'svg', True),
            #     "data": {},
            #     "children": []
            # })
            # _products_dict['data'].update(
            #     self.popupBuilder({"smiles": _parent, "generation": "0"},
            #         self.metabolite_keys, "{}".format(self.metID),
            #         "Metabolite Information"
            #     )
            # )
            # self.products_list.append(_parent)

            # skipping 2nd parent metabolite:
            second_parent = root[_parent]['metabolites'].keys()[0]
            root = root[_parent]['metabolites'][second_parent]
            # not-skipping version without 2nd parent problem:
            # root = root[_parent]

            _products_dict.update({
                "id": self.metID,
                "name": "<img class='blank_node' src='/static_qed/cts/images/loader_node.gif' />",
                # "name": self.nodeWrapper(_parent, self.tree_image_height, self.tree_image_width, self.image_scale, self.metID,'svg', True),
                # 'name': "",
                "data": {'smiles': _parent, 'routes': root['routes'], 'generation': root['generation']},
                "children": []
            })

            
        else:
            if root['generation'] > 0 and root['generation'] <= gen_limit:
                # continue walking tree until generation limit is met..
                # _products_dict.update({
                #     "id": self.metID,
                #     "name": self.nodeWrapper(root['smiles'], self.tree_image_height, self.tree_image_width, self.image_scale, self.metID, 'svg', True),
                #     "data": {},
                #     "children": []
                # })
                # _products_dict['data'].update(self.popupBuilder(root, self.metabolite_keys, "{}".format(self.metID), "Metabolite Information"))
                _products_dict.update({
                    "id": self.metID,
                    "name": "<img class='blank_node' src='/static_qed/cts/images/loader_node.gif' />",
                    # 'name': self.nodeWrapper(root['smiles'], self.tree_image_height, self.tree_image_width, self.image_scale, self.metID, 'svg', True),
                    # 'name': "",
                    "data": {'smiles': root['smiles'], 'routes': root['routes'].split(',')[-1], 'generation': root['generation']},
                    "children": []
                })
                # self.products_list.append(root['smiles'])


        for key, value in root.items():
            if isinstance(value, dict):
                for key2, value2 in root[key].items():
                    root2 = root[key][key2]
                    # if len(root2) > 0 and 'children' in _products_dict and root2['generation'] < gen_limit:
                    if len(root2) > 0 and 'children' in _products_dict and root['generation'] < gen_limit:
                        # continue walking branch if root2 has contents, one of those contents is 'children', and
                        # the generation limit isn't exceeded..
                        _products_dict['children'].append(self.traverse(root2, gen_limit))

        return _products_dict



    def make_data_request(self, structure, prop_obj, method=None):
        # url = self.baseUrl + prop_obj.url
        url = self.efs_server_url + self.efs_metabolizer_endpoint

        _valid_result = False  # for retry logic
        _retries = 0
        while not _valid_result and _retries < self.max_retries:
            # retry data request to chemaxon server until max retries or a valid result is returned
            try:
                response = requests.post(url, data=json.dumps(post_data), headers=self.headers, timeout=self.request_timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    prop_obj.results = json.loads(response.content)
                    return prop_obj.results
                _retries += 1
            except Exception as e:
                logging.warning("Exception in metabolizer_calculator.py: {}".format(e))
                _retries += 1

            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))


    def data_request_handler(self, request_dict):

        logging.warning("$$$ METABOLIZER REQUEST: {} $$$".format(request_dict))

        # reactionLibs = {
        #     "hydrolysis": request_dict.get('abiotic_hydrolysis'),
        #     "abiotic_reduction": request_dict.get('abiotic_reduction'),
        #     # "human_biotransformation": self.mamm_metabolism
        # }

        # _trans_libs = []
        # for key, value in reactionLibs.items():
        #     if value:
        #         _trans_libs.append(key)

        # if not request_dict.get('gen_limit', False):
        #     request_dict['gen_limit'] = 1

        # # NOTE: populationLimit is hard-coded to 0 as it currently does nothing

        # _data_dict = {
        #     'structure': request_dict['chemical'],
        #     'generationLimit': request_dict['gen_limit'],
        #     'populationLimit': 0,
        #     # 'likelyLimit': self.likely_limit,
        #     'likelyLimit': 0.001,
        #     'excludeCondition': ""  # 'generateImages': False
        # }

        # if len(_trans_libs) > 0:
        #     _data_dict.update({'transformationLibraries': _trans_libs})

        _data_dict = request_dict.get('metabolizer_post')
        logging.info("METABOLIZER POST: {}".format(_data_dict))

        response = self.getTransProducts(_data_dict)
        # response = self.make_data_request(request_dict['chemical'], self, None)
        _results = MetabolizerCalc().recursive(response, int(request_dict['gen_limit']))

        _products_data = json.loads(_results)

        _response_obj = {
            'calc': "chemaxon",  # todo: change to metabolizer, change in template too
            'prop': "products",
            'node': request_dict.get('node'),
            # 'data': json.loads(_results),
            'data': _products_data['tree'],
            'total_products': _products_data['total_products'],
            'chemical': request_dict.get('chemical'),
            'workflow': 'gentrans',
            'run_type': 'batch',
            'request_post': request_dict            
        }


        return _response_obj


    def validate_response(self, response):
        """
        Validates jchem response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("cts_celery calculator_chemaxon -- jchem server response status not 200, but: {}".format(response.status_code))
            # raise Exception("non 200 response from jchem: {}".format(response))
            return False
        
        # successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
        # json_obj = json.loads(response.content)

        # TODO: verify if blank data, finding the source of the empty water sol values...
        
        return True


    def getTransProducts(self, request_obj):
        """
        Makes request to metabolizer
        """
        url = self.efs_server_url + self.efs_metabolizer_endpoint
        self.request_timeout = 60
        return self.web_call(url, request_obj)