import requests
import json
import logging
import os
# import redis
from .calculator import Calculator



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

        self.transformation_libraries = [
            "hydrolysis",
            "abiotic_reduction",
            "photolysis_unranked",
            "photolysis_ranked",
            "mammalian_metabolism",
            "combined_abioticreduction_hydrolysis",
            "combined_photolysis_abiotic_hydrolysis"
        ]

        # CTSWS Transformation Products Request
        self.transformation_request = {
            'structure': None,
            'generationLimit': 1,  # make sure to get this from front end
            'populationLimit': 0,
            'likelyLimit': 0.1,
            'transformationLibraries': ["hydrolysis", "abiotic_reduction"],  # NOTE: no transformationLibraries key:val for mammalian metabolism
            'excludeCondition': "hasValenceError()"
        }

        self.response_obj = {
            'calc': "metabolizer",  # todo: change to metabolizer, change in template too
            'prop': "products",
            'node': None,
            'data': None,
            'total_products': 0,
            'unique_products': 0,
            'chemical': None,
            'workflow': "gentrans",
            'run_type': None,
            'request_post': None
        }

        self.unique_products = []


    def recursive(self, jsonDict, gen_limit, unranked=False):
        """
        Starting point for walking through
        metabolites dictionary and building json
        that thejit (visualization javascript
        library) understands
        """
        root = jsonDict['results']

        reDict = {}
        reDict.update({
            'tree': self.traverse(root, gen_limit, unranked),
            'total_products': self.metID - 1,  # subtract out the parent for "total products" value
            'unique_products': len(self.unique_products)
        })

        self.metID = 0  # resets the metID attribute
        self.unique_products = []

        # Need to hit SMILES filter, then retrieve molecular info,
        # node image, and popup image for products...
        # logging.info("PRODUCTS LIST: {}".format(self.products_list))

        
        # Could make one big request to get node images from list of 
        # nodes here! But, the tree would need to be walked again so it
        # can have edit the 'name' key with the node images...


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

            _parent = list(root.keys())[0]  # python 3 fix
        
            # root = root[_parent]['metabolites']  # ctsws update
            root = root[_parent][0]  # ctsws update 2
            
            _products_dict.update({
                "id": self.metID,
                "name": "<img class='blank_node' src='" + os.getenv("CTS_STATIC_URL") + "cts_app/images/loader_node.gif' />",
                "data": {
                    'smiles': _parent,
                    'routes': root['route'],
                    'generation': root['generation'],
                    # 'accumulation': round(root.get('accumulation'), 4),
                    'accumulation': "N/A" if unranked else round(100 * root.get('accumulation'), 2),
                    'production': "N/A" if unranked else round(100 * root.get('production'), 2),
                    'globalAccumulation': "N/A" if unranked else round(100 * root.get('globalAccumulation'), 2),
                    'likelihood': "N/A" if unranked else root.get('likelihood')
                },
                "children": []
            })
            
        else:

            if root['generation'] > 0 and root['generation'] <= gen_limit:

                likelihood = self.setLikelyhoodValue(root)

                # continue walking tree until generation limit is met..
                _products_dict.update({
                    "id": self.metID,
                    "name": "<img class='blank_node' src='" + os.getenv("CTS_STATIC_URL") + "cts_app/images/loader_node.gif' />",
                    "data": {
                        'smiles': root['smiles'],
                        'routes': root['route'].split(',')[-1],
                        'generation': root['generation'],
                        'accumulation': "N/A" if unranked else round(100 * root.get('accumulation'), 2),
                        'production': "N/A" if unranked else round(100 * root.get('production'), 2),
                        'globalAccumulation': "N/A" if unranked else round(100 * root.get('globalAccumulation'), 2),
                        'likelihood': "N/A" if unranked else likelihood
                    },
                    "children": []
                })

                # UNIQUE PRODUCTS DETERMINED HERE
                # self.products_list.append(root['smiles'])
                if not root['smiles'] in self.unique_products:
                    self.unique_products.append(root['smiles'])

        # for key, value in root.items():
        #     if isinstance(value, list) and len(value) > 0:
        #         for met_obj in value:
        #             if 'children' in _products_dict and root['generation'] < gen_limit:
        #                 _products_dict['children'].append(self.traverse(met_obj, gen_limit, unranked))

        # for key, value in root.items():
        # if isinstance(value, list) and len(value) > 0:
        if 'metabolites' in root and len(root['metabolites']) > 0:
            for met_obj in root['metabolites']:
                if 'children' in _products_dict and root['generation'] < gen_limit:
                    _products_dict['children'].append(self.traverse(met_obj, gen_limit, unranked))

        return _products_dict



    def data_request_handler(self, request_dict):

        _data_dict = request_dict.get('metabolizer_post')
        _data_dict.update({'structure': request_dict.get('chemical'), 'excludeCondition': 'hasValenceError()'})

        _response_obj = dict(self.response_obj)
        _response_obj['node'] = request_dict.get('node'),
        _response_obj['chemical'] = request_dict.get('chemical'),
        _response_obj['run_type'] = request_dict.get('run_type', 'batch'),
        _response_obj['request_post'] = request_dict

        unranked = False
        if 'photolysis_unranked' in _data_dict.get('transformationLibraries', []):
            unranked = True

        response = self.getTransProducts(_data_dict)

        if "error" in response:
            _response_obj["error"] = response["error"]
            return _response_obj
        elif response.get("status") == "error":
            _response_obj["error"] = "Error getting transformation products"
            return _response_obj

        _results = self.recursive(response, int(request_dict['gen_limit']), unranked)

        _products_data = json.loads(_results)

        _response_obj['data'] = _products_data['tree']
        _response_obj['total_products'] = _products_data['total_products']
        _response_obj['unique_products'] = _products_data['unique_products']

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
        structure = request_obj.get("structure")
        gen_limit = request_obj.get("generationLimit")
        trans_libs = request_obj.get("transformationLibraries", [])

        if gen_limit > 4:
            _response_obj = dict(request_obj)
            _response_obj["error"] = "Must request generation limit <= 4 generations."
            return _response_obj

        if len(trans_libs) != 1:
            _response_obj = dict(request_obj)
            _response_obj["error"] = "Can only run one transformation library at a time."
            return _response_obj

        url = self.efs_server_url + self.efs_metabolizer_endpoint
        self.request_timeout = 120
        return self.web_call(url, request_obj)



    def setLikelyhoodValue(self, product_data):
        """
        Checks likelihood value from CTSWS, keeps value as
        "LIKELY" if likelihood > 10%, "UNLIKELY" if < 0.1%, and
        sets likelihood as "PROBABLE" if it's between 0.1% and 10%.
        """
        global_accumulation = product_data.get('globalAccumulation')

        if global_accumulation < 0.001:
            return "UNLIKELY"
        elif global_accumulation > 0.001 and global_accumulation < 0.1:
            return "PROBABLE"
        elif global_accumulation > 0.1:
            return "LIKELY"
        else:
            return product_data.get('likelihood')