import requests
import json
import logging
import os
import redis

import smilesfilter
# import jchem_properties
# import jchem_rest
# import data_walks
from calculator import Calculator


class MetabolizerCalc(Calculator):
    """
    Calc class for CTSWS (formerly EFS) Metabolizer.

    This was originally handled by ChemaxonCalc and data_walks,
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

        # CTSWS Transformation Products Request
        self.transformation_request = {
            'structure': None,
            'generationLimit': 1,  # make sure to get this from front end
            'populationLimit': 0,
            'likelyLimit': 0.001,
            'transformationLibraries': ["hydrolysis", "abiotic_reduction"],  # NOTE: no transformationLibraries key:val for mammalian metabolism
            'excludeCondition': ""
        }


    # def data_request_handler(self, request_dict, message=None):
    #     """
    #     This will basically be data_walks main code.
    #     """




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
        reDict.update(traverse(root, gen_limit))
        return json.dumps(reDict)


    metID = 0  # unique id for each node
    metabolite_keys = ['smiles', 'formula', 'iupac', 'mass', 'accumulation', 'production', 'transmissivity', 'generation', 'routes', 'exactMass']
    image_scale = 50

    def traverse(self, root, gen_limit):
        """
        For gentrans model output - metabolites tree
        """

        global metID
        metID += 1
        newDict = {}

        logging.info("metabolites: {}".format(metID))

        tblID = "{}_table".format(metID)  # id for node's tooltip table

        if metID == 1:
            parent = root.keys()[0]
            newDict.update({"id": metID, "name": nodeWrapper(parent, 114, 100, image_scale, metID, 'svg', True), "data": {}, "children": []})
            # newDict.update({"id": metID, "name": nodeWrapper(parent, None, 100, 28), "data": {}, "children": []})
            newDict['data'].update(popupBuilder({"smiles": parent, "generation": "0"}, metabolite_keys, "{}".format(metID),
                                                "Metabolite Information"))

            # request_obj = {'chemical': parent}  # chemical info request object
            # mol_info = jchem_rest.getChemDetails(request_obj)
            filtered_smiles = smilesfilter.filterSMILES(parent)
            mol_info = jchem_rest.getChemDetails({'chemical': filtered_smiles})
            
            if 'data' in mol_info:
                for key, val in mol_info['data'][0].items():
                    if key in metabolite_keys:
                        newDict['data'].update({key: val})

            # skipping 2nd parent metabolite:
            second_parent = root[parent]['metabolites'].keys()[0]
            root = root[parent]['metabolites'][second_parent]

            # not-skipping version without 2nd parent problem:
            # root = root[parent]
            
        else:
            if root['generation'] > 0 and root['generation'] <= gen_limit:
                newDict.update({"id": metID, "name": nodeWrapper(root['smiles'], 114, 100, image_scale, metID, 'svg', True), "data": {}, "children": []})
                # newDict.update({"id": metID, "name": nodeWrapper(root['smiles'], None, 100, 28), "data": {}, "children": []})
                newDict['data'].update(popupBuilder(root, metabolite_keys, "{}".format(metID), "Metabolite Information"))

                # request_obj = {'chemical': root['smiles']}
                # mol_info = jchem_rest.getChemDetails(request_obj)
                filtered_smiles = smilesfilter.filterSMILES(root['smiles'])
                mol_info = jchem_rest.getChemDetails({'chemical': filtered_smiles})
                
                if 'data' in mol_info:
                    for key, val in mol_info['data'][0].items():
                        if key in metabolite_keys:
                            newDict['data'].update({key: val})

        for key, value in root.items():
            if isinstance(value, dict):
                for key2, value2 in root[key].items():
                    root2 = root[key][key2]
                    if len(root2) > 0 and 'children' in newDict and root['generation'] < gen_limit:
                    # if len(root2) > 0 and 'children' in newDict and root2['generation'] < gen_limit:
                        newDict['children'].append(traverse(root2, gen_limit))

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

        logging.info("JCHEM REQUEST URL: {}".format(url))
        logging.info("JCHEM REQUEST POST: {}".format(post_data))

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
                logging.warning("Exception in jchem_calculator.py: {}".format(e))
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