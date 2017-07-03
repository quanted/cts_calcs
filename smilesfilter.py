__author__ = 'KWOLFE'

import requests
import logging
import json
# import jchem_rest
from .calculator import Calculator
from .jchem_properties import Tautomerization
# from calculator_chemaxon import JchemProperty
# from calculator_chemaxon import Tautomerization

max_weight = 1500 # max weight [g/mol] for epi, test, and sparc


def is_valid_smiles(smiles):

    excludestring = {".","[Ag]","[Al]","[Au]","[As]","[As+","[B]","[B-]","[Br-]","[Ca]",
                        "[Ca+","[Cl-]","[Co]","[Co+","[Fe]","[Fe+","[Hg]","[K]","[K+","[Li]",
                        "[Li+","[Mg]","[Mg+","[Na]","[Na+","[Pb]","[Pb2+]","[Pb+","[Pt]",
                        "[Sc]","[Si]","[Si+","[SiH]","[Sn]","[W]"}

    return_val = {
        "valid" : False,
        "smiles": smiles,
        "processedsmiles" : ""
    }

    if any(x in smiles for x in excludestring):
        return return_val

    try:
        processed_smiles = filterSMILES(smiles)
    except Exception as e:
        logging.warning("!!! Error in smilesfilter {} !!!".format(e))
        raise "smiles filter exception, possibly invalid smiles..."
            
    return_val["valid"] = True
    return_val["processedsmiles"] = processed_smiles

    return return_val


# def filterSMILES(smiles):
#     """
#     calculator-independent SMILES processing.
#     uses jchem web services through jchem_rest
#     """
#     response = jchem_rest.filterSMILES({'smiles': smiles})
#     logging.warning("FILTER RESPONSE: {}".format(response))
#     try:
#         filtered_smiles = response['results'][-1] # picks out smiles from efs???
#         logging.warning("NEW SMILES: {}".format(filtered_smiles))
#         return filtered_smiles
#     except Exception as e:
#         logging.warning("> error in filterSMILES: {}".format(e))
#         raise e
def singleFilter( request_obj):
    """
    Calls single EFS Standardizer filter
    for filtering SMILES
    """
    try:
        smiles = request_obj.get('smiles')
        action = request_obj.get('action')
    except Exception as e:
        logging.info("Exception retrieving mass from jchem: {}".format(e))
        raise
    post_data = {
        "structure": smiles,
        "actions": [
            action
        ]
    }
    url = Calculator().efs_server_url + Calculator().efs_standardizer_endpoint
    return Calculator().web_call(url, post_data)


# def filterSMILES(request_obj):
def filterSMILES(smiles):
    """
    cts ws call to jchem to perform various
    smiles processing before being sent to
    p-chem calculators
    """

    # smiles = request_obj.get('smiles')

    # Updated approach (todo: more efficient to have CTSWS use major taut instead of canonical)
    # 1. CTSWS actions "removeExplicitH" and "transform".
    url = Calculator().efs_server_url + Calculator().efs_standardizer_endpoint
    post_data = {
        'structure': smiles,
        'actions': [
            "removeExplicitH",
            "transform"
        ]
    }
    response = Calculator().web_call(url, post_data)

    filtered_smiles = response['results'][-1] # picks last item, format: [filter1 smiles, filter1 + filter2 smiles]
    
    # 2. Get major tautomer from jchem:
    # tautObj = JchemCalc.getPropObject('tautomerization')
    # tautObj.setPostDataValues({"calculationType": "MAJOR"})
    taut_obj = Tautomerization()
    taut_obj.postData.update({'calculationType': 'MAJOR'})
    # taut_obj.makeDataRequest(filtered_smiles)
    taut_obj.make_data_request(filtered_smiles, taut_obj)

    # try:
    #     _taut_response = Calculator().web_call(taut_obj.url, taut_obj.postData)
    #     major_taut_smiles = _taut_response['result']['structureData']['structure']
    # except Exception as e:
    #     logging.warning("filterSMILES excpetion: {}".format(e))


    # todo: verify this is major taut result smiles, not original smiles for major taut request...
    major_taut_smiles = taut_obj.results['result']['structureData']['structure']

    # logging.warning("MAJOR TAUT SMILES: {}".format(major_taut_smiles))

    # 3. Using major taut smiles for final "neutralize" filter:
    post_data = {
        'structure': major_taut_smiles, 
        'actions': [
            "neutralize"
        ]
    }
    response = Calculator().web_call(url, post_data)

    final_smiles = response['results'][-1]
    logging.warning("FINAL FITERED SMILES: {}".format(final_smiles))

    return response


def checkMass(chemical):
    """
    returns true if chemical mass is less
    than 1500 g/mol
    """
    logging.info("checking mass..")
    try:
        json_obj = Calculator().getMass({'chemical': chemical}) # get mass from jchem ws
    except Exception as e:
        logging.warning("!!! Error in checkMass() {} !!!".format(e))
        raise e
    logging.info("mass response data: {}".format(json_obj))
    struct_mass = json_obj['data'][0]['mass']
    logging.info("structure's mass: {}".format(struct_mass))

    if struct_mass < 1500  and struct_mass > 0:
        return True
    else:
        return False


def clearStereos(smiles):
    """
    clears stereoisomers from smiles
    """
    try:
        response = singleFilter({'smiles':smiles, 'action': "clearStereo"})
        filtered_smiles = response['results'] # get stereoless structure
    except Exception as e:
        logging.warning("!!! Error in clearStereos() {} !!!".format(e))
        raise e
    return filtered_smiles


def transformSMILES(smiles):
    """
    N(=O)=O >> [N+](=O)[O-]
    """
    try:
        response = singleFilter({'smiles':smiles, 'action': "transform"})
        filtered_smiles = response['results'] # get stereoless structure
    except Exception as e:
        logging.warning("!!! Error in transformSMILES() {} !!!".format(e))
        raise e
    return filtered_smiles


def untransformSMILES(smiles):
    """
    [N+](=O)[O-] >> N(=O)=O
    """
    try:
        response = singleFilter({'smiles':smiles, 'action': "untransform"})
        filtered_smiles = response['results'] # get stereoless structure
    except Exception as e:
        logging.warning("!!! Error in untransformSMILES() {} !!!".format(e))
        raise e
    return filtered_smiles


def parseSmilesByCalculator(structure, calculator):
    """
    Calculator-dependent SMILES filtering!
    """

    logging.info("Parsing SMILES by calculator..")
    filtered_smiles = structure

    #1. check structure mass..
    if calculator != 'chemaxon':
        logging.info("checking mass for: {}...".format(structure))
        if not checkMass(structure):
            logging.info("Structure too large, must be < 1500 g/mol..")
            raise "Structure too large, must be < 1500 g/mol.."

    #2-3. clear stereos from structure, untransform [N+](=O)[O-] >> N(=O)=O..
    if calculator == 'epi' or calculator == 'sparc' or calculator == 'measured':
        try:
            # clear stereoisomers:
            filtered_smiles = clearStereos(structure)
            logging.info("stereos cleared: {}".format(filtered_smiles))

            # transform structure:
            filtered_smiles = str(filtered_smiles[-1])
            filtered_smiles = str(untransformSMILES(filtered_smiles)[-1])
            logging.info("structure transformed..")
        except Exception as e:
            logging.warning("!!! Error in parseSmilesByCalculator() {} !!!".format(e))
            raise e

    # 4. Check for metals and stuff (square brackets):
    if calculator == 'epi' or calculator == 'measured':
        if '[' in filtered_smiles or ']' in filtered_smiles:
            # bubble up to calc for handling error
            raise Exception("{} cannot process metals...".format(calculator))
            # logging.warning("EPI ignoring request due to brackets in SMILES..")
            # postData.update({'data': "EPI Suite cannot process charged species or metals (e.g., [S+], [c+])"})
            # if redis_conn and sessionid:
            #     for prop in props:
            #         postData['prop'] = prop
            #         postData['node'] = node
            #         if run_type:
            #             postData['run_type'] = run_type
            #         redis_conn.publish(sessionid, json.dumps(postData))
            #     return

    return filtered_smiles