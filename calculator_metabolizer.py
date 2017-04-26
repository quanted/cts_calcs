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
        parent = root.keys()[0]
        reDict.update(self.traverse(root, gen_limit))
        return json.dumps(reDict)


    def traverse(self, root, gen_limit):
        """
        For gentrans model output - metabolites tree
        """

        # global metID
        self.metID += 1
        newDict = {}

        logging.info("metabolites: {}".format(self.metID))

        tblID = "{}_table".format(self.metID)  # id for node's tooltip table

        if self.metID == 1:
            parent = root.keys()[0]
            newDict.update({"id": self.metID, "name": self.nodeWrapper(parent, self.tree_image_height, self.tree_image_width, self.image_scale, self.metID, 'svg', True), "data": {}, "children": []})
            newDict['data'].update(self.popupBuilder({"smiles": parent, "generation": "0"}, self.metabolite_keys, "{}".format(self.metID),
                                                "Metabolite Information"))

            _filtered_smiles = smilesfilter.filterSMILES(parent)['results'][-1]
            _mol_info = self.getChemDetails({'chemical': _filtered_smiles})
            
            if 'data' in _mol_info:
                for key, val in _mol_info['data'][0].items():
                    if key in self.metabolite_keys:
                        newDict['data'].update({key: val})

            # skipping 2nd parent metabolite:
            second_parent = root[parent]['metabolites'].keys()[0]
            root = root[parent]['metabolites'][second_parent]
            # not-skipping version without 2nd parent problem:
            # root = root[parent]
            
        else:
            if root['generation'] > 0 and root['generation'] <= gen_limit:
                newDict.update({"id": self.metID, "name": self.nodeWrapper(root['smiles'], self.tree_image_height, self.tree_image_width, self.image_scale, self.metID, 'svg', True), "data": {}, "children": []})
                newDict['data'].update(self.popupBuilder(root, self.metabolite_keys, "{}".format(self.metID), "Metabolite Information"))

                _filtered_smiles = smilesfilter.filterSMILES(root['smiles'])['results'][-1]
                _mol_info = self.getChemDetails({'chemical': _filtered_smiles})
                
                if 'data' in _mol_info:
                    for key, val in _mol_info['data'][0].items():
                        if key in self.metabolite_keys:
                            newDict['data'].update({key: val})

        for key, value in root.items():
            if isinstance(value, dict):
                for key2, value2 in root[key].items():
                    root2 = root[key][key2]
                    if len(root2) > 0 and 'children' in newDict and root['generation'] < gen_limit:
                    # if len(root2) > 0 and 'children' in newDict and root2['generation'] < gen_limit:
                        newDict['children'].append(self.traverse(root2, gen_limit))

        return newDict



    def make_data_request(self, structure, prop_obj, method=None):
        url = self.baseUrl + prop_obj.url
        prop_obj.postData.update({
            "result-display": {
                "include": ["structureData", "image"],
                "parameters": {
                    "structureData": "smiles"
                }
            }
        })
        post_data = {
            "structure": structure,
            "parameters": prop_obj.postData
        }

        logging.info("EFS REQUEST URL: {}".format(url))
        logging.info("EFS REQUEST POST: {}".format(post_data))

        if method:
            post_data['parameters']['method'] = method


        _valid_result = False  # for retry logic
        _retries = 0
        while not _valid_result and _retries < self.max_retries:
            # retry data request to chemaxon server until max retries or a valid result is returned
            try:
                response = requests.post(url, data=json.dumps(post_data), headers=self.headers, timeout=self.request_timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    # self.results = json.loads(response.content)
                    prop_obj.results = json.loads(response.content)
                    # logging.info("Response from jchem server: {}".format(prop_obj.results))
                    break
                _retries += 1
            except Exception as e:
                logging.warning("Exception in metabolizer_calculator.py: {}".format(e))
                _retries += 1

            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))


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
        return self.web_call(url, request_obj)