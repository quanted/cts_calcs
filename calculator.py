__author__ = 'np'

import requests
import json
import logging
import os
import redis

import smilesfilter
# from chemaxon_cts import jchem_properties
# from chemaxon_cts import jchem_calculator
# from epi_cts import epi_calculator
# from sparc_cts import sparc_calculator
# import calculator_chemaxon


headers = {'Content-Type': 'application/json'}

class Calculator(object):
    """
    Skeleton class for calculators
    """

    def __init__(self, calc=None):
        """
        calc -- p-chem calculator name, 'getTransProducts', or 'getSpeciationData'
        """
        self.name = ''
        self.propMap = {}
        self.baseUrl = None
        self.urlStruct = ''
        self.results = ''
        self.headers = {'Content-Type': 'application/json'}
        self.request_timeout = 30  # default, set unique ones in calc sub classes
        self.max_retries = 3
        # self.request_type = ''  # http or ws (https or wss [todo])

        self.redis_hostname = os.environ.get('REDIS_HOSTNAME')
        self.redis_port = os.environ.get('REDIS_PORT')
        self.redis_conn = redis.StrictRedis(host=self.redis_hostname, port=self.redis_port, db=0)

        # cts p-chem properties
        self.pchem_props = [
            'boiling_point',
            'melting_point',
            'water_sol',
            'vapor_press',
            'mol_diss',
            'ion_con',
            'henrys_law_con',
            'kow_no_ph',
            'kow_wph',
            'kow_ph',
            'kow'
        ]

        # cts chemical information dict
        self.chemical_information = {
            'chemical': None,  # user-entered chemical (as-entered or drawn)
            'orig_smiles': None,  # original conversion to SMILES
            'smiles': None,  # SMILES after filtering, used for calculations
            'formula': None,
            'iupac': None,
            'mass': None,
            'structureData': None,  # drawn chemical structure format for MarvinSketch
            'exactMass': None,
        }

        # cts chemical information request
        self.chemical_information_request = {
            'chemical': None,
            'get_structure_data': False,
        }

        # cts api data object for p-chem data request
        self.data_obj = {
            'calc': None,
            'prop': None,
            'data': None,
            'chemical': None,
        }

        # cts p-chem request object with default key:vals.
        # can handle list of props (ws) or single prop (cts api)
        self.pchem_request = {
            'service': None,
            'chemical': None,
            'prop': None,
            'sessionid': None,
            'method': None,
            'ph': 7.0,
            'node': None,
            'calc': None,
            'run_type': None,
            'workflow': None,
            'mass': None,
            'props': [],
        }

        # cts p-chem response object with defaults, returns one prop per reponse
        self.pchem_response = {
            'chemical': None,
            'calc': None,
            'prop': None,
            'method': None,
            'run_type': None,
            'workflow': None,
            'node': None,
            'request_post': None,
            'data': None,
            'error': False,  # ehh?
        }

        if calc:
            self.getCalcObject(calc)


    @classmethod
    def getCalcObject(self, calc):
        """
        Returns instance of calculator:
        chemaxon, epi, test, sparc, or measured
        """
        if calc == 'chemaxon':
            # return calculator_chemaxon.ChemaxonCalc()
            return None
        # elif calc == 'epi':
            # return EpiCalc()
        # elif calc == 'sparc':
        #     return sparc_calculator.SparcCalc()
        else:
            return None


    def getUrl(self, prop):
        if prop in self.propMap:
            calcProp = self.propMap[prop]['urlKey']
            return self.urlStruct.format(calcProp)
        else:
            return "Error: url key not found"

    def getPropKey(self, prop):
        if prop in self.propMap:
            return self.propMap[prop]['propKey']
        else:
            return "Error: prop not found"

    def getResultKey(self, prop):
        if prop in self.propMap:
            return self.propMap[prop]['resultKey']
        else:
            return "Error: result key not found"