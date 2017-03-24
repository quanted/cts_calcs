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

        # chemical information dict
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
            'props': None,
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

    # def request_logic(self, request, message=None): 
    #     """
    #     request_manager generalized.

    #     Takes in request dict and makes requests
    #     to calculator servers. Also filters smiles based
    #     on calc being requested, and handles exceptions when
    #     making requests to calc servers.

    #     message - django channels message object, if there is one
    #     then it's WS protocol and reply_channel for data push

    #     Returns response dict for ws/consumers or
    #     cts_api/http to consume.

    #     Note: Takes errors raised from calc classes and sends
    #     them as 'data' for error response to frontend.
    #     Note: errors could also raise from here to consumers
    #     or cts_api endpoints if that makes more sense.
    #     """

    #     _request_dict = self.pchem_request  # request dict
    #     _response = self.pchem_response  # response dict

    #     # Fill calculator p-chem request dict with incoming request values
    #     for key in _request_dict.keys():
    #         _request_dict[key] = request.get(key)


    #     # Throughout the years, these are various ways that getting
    #     # the 'props' list was performed.
    #     # TODO: Move this to cts_rest? Or wherever the incoming request is??
    #     # try:
    #     #     props = request.POST.get("props[]")
    #     #     if not props:
    #     #         props = request.POST.getlist("props")
    #     # except AttributeError:
    #     #     props = request.POST.get("props")


    #     # Filtering chemical (todo: keep track of chemical, orig smiles, filtered smiles)
    #     filtered_smiles = ''
    #     try:
    #         filtered_smiles = smilesfilter.parseSmilesByCalculator(_request_dict['chemical'], _request_dict['calc']) # call smilesfilter
    #         # _request_dict['chemical'] = smilesfilter.parseSmilesByCalculator(_request_dict['chemical'], _request_dict['calc']) # call smilesfilter
    #     except Exception as err:
    #         logging.warning("Error filtering SMILES: {}".format(err))
    #         _request_dict.update({'data': 'Cannot filter SMILES for EPI data'})
    #         return _request_dict


    #     # get calc object based on calc's name:
    #     calc_obj = self.getCalcObject(_request_dict['calc'])  # calculator.py class func

    #     # todo: remove this as requirement for api docs calls
    #     if _request_dict['run_type'] == 'rest':
    #         _request_dict['props'] = [_request_dict['prop']]  # rest api currently does single prop calls

    #     logging.info("Request received by calculator request_logic: {}".format(_request_dict))

    #     # call calc subclass make_data_request(), which should take care of the
    #     # logic below as it varies across calcs.
    #     _response_dict = {}
    #     try:
    #         _response_dict = calc_obj.data_request_handler(_request_dict, message)  # call calc-specific data request handler
    #     except Exception as e:
    #         logging.warning("Exception in calculator.py: {}".format(e))
    #         _response_dict['data'] = e  # assuming calc sub classes are sending custom excpetions for user..

    #     # return response of data for single prop
    #     return _response_dict





class SparcCalc(Calculator):
    def __init__(self, smiles=None, meltingpoint=0.0, pressure=760.0, temperature=25.0):

        self.base_url = os.environ['CTS_SPARC_SERVER']
        self.multiproperty_url = '/sparc-integration/rest/calc/multiProperty'
        self.name = "sparc"
        self.smiles = smiles
        self.solvents = dict()
        self.pressure = pressure
        self.meltingpoint = meltingpoint
        self.temperature = temperature
        self.propMap = {
            "water_sol" : "SOLUBILITY",
            "vapor_press" : "VAPOR_PRESSURE",
            "henrys_law_con" : "HENRYS_CONSTANT",
            "mol_diss" : "DIFFUSION",
            "boiling_point": "BOILING_POINT"
        }
        self.sparc_props = {
            "SOLUBILITY": "water_sol",
            "VAPOR_PRESSURE": "vapor_press",
            "WATER_DIFFUSION": "mol_diss",
            "FULL_SPECIATION": "ion_con",
            "HENRYS_CONSTANT": "henrys_law_con",
            "DISTRIBUTION": "kow_no_ph",
            "LOGD": "kow_wph",
            "BOILING_POINT": "boiling_point" 
        }

    def get_sparc_query(self):
        query = {
            'pressure': self.pressure,
            'meltingPoint': self.meltingpoint,
            'temperature': self.temperature,
            'calculations': self.getCalculations(),
            'smiles': self.smiles,
            'userId': None,
            'apiKey': None,
            'type': 'MULTIPLE_PROPERTY',
            'doSolventInit': False
        }
        return query

    def get_calculation(self, sparc_prop=None, units=None):
        calc = {
            'solvents': [],
            'units': units,
            'pressure': self.pressure,
            'meltingPoint': self.meltingpoint,
            'temperature': self.temperature,
            'type': sparc_prop
        }
        return calc

    def get_solvent(self, smiles=None, name=None):
        solvent = {
            'solvents': None,
            'smiles': smiles,
            'mixedSolvent': False,
            'name': name
        }
        return solvent


    def getCalculations(self):
        calculations = []
        calculations.append(self.get_calculation("VAPOR_PRESSURE", "Torr"))
        calculations.append(self.get_calculation("BOILING_POINT", "degreesC"))
        calculations.append(self.get_calculation("DIFFUSION", "NO_UNITS"))
        calculations.append(self.get_calculation("VOLUME", "cmCubedPerMole"))
        calculations.append(self.get_calculation("DENSITY", "gPercmCubed"))
        calculations.append(self.get_calculation("POLARIZABLITY", "angCubedPerMolecule"))
        calculations.append(self.get_calculation("INDEX_OF_REFRACTION", "dummy"))

        calcHC = self.get_calculation("HENRYS_CONSTANT", "AtmPerMolPerM3")
        calcHC["solvents"].append(self.get_solvent("O", "water"))
        calculations.append(calcHC)

        calcSol = self.get_calculation("SOLUBILITY", "mgPerL")
        calcSol["solvents"].append(self.get_solvent("O", "water"))
        calculations.append(calcSol)

        calcAct = self.get_calculation("ACTIVITY", "dummy")
        calcAct["solvents"].append(self.get_solvent("O", "water"))
        calculations.append(calcAct)

        calculations.append(self.get_calculation("ELECTRON_AFFINITY", "dummy"))

        calcDist = self.get_calculation("DISTRIBUTION", "NO_UNITS")
        calcDist["solvents"].append(self.get_solvent("O", "water"))
        calcDist["solvents"].append(self.get_solvent("OCCCCCCCC", "octanol"))

        calculations.append(calcDist)

        return calculations


    def makeDataRequest(self):

        post = self.get_sparc_query()
        url = self.base_url + self.multiproperty_url

        logging.info("SPARC URL: {}".format(url))
        logging.info("SPARC POST: {}".format(post))

        try:
            response = requests.post(url, data=json.dumps(post), headers=headers, timeout=20)
            self.results = json.loads(response.content)
        except requests.exceptions.ConnectionError as ce:
            logging.info("connection exception: {}".format(ce))
            raise
        except requests.exceptions.Timeout as te:
            logging.info("timeout exception: {}".format(te))
            raise
        else:
            return self.results


    def parseMultiPropResponse(self, results):
        """
        Loops through data grabbing the results
        and building {calc, prop, data} objects
        for front end.

        TODO: Add more info to returned data (object instead of value)
        """
        if not results or not isinstance(results, list):
            raise Exception("(sparc) no results given to parse or not a list..")

        sparc_response = []
        for prop_data in results:
            sparc_prop = prop_data['type']
            if sparc_prop in self.sparc_props:
                cts_prop = self.sparc_props[sparc_prop]
                data = prop_data['result']
                sparc_response.append({'calc': 'sparc', 'prop': cts_prop, 'data': data})

        return sparc_response


    def makeCallForPka(self):
        """
        Separate call for pKa. Not sure why but
        what what I'm told it needs to be done 
        separately for now
        """
        pka_url = "/sparc-integration/rest/calc/fullSpeciation"
        url = self.base_url + pka_url
        logging.info("URL: {}".format(url))
        sparc_post = {
            "type":"FULL_SPECIATION",
            "temperature":25.0,
            "minPh":0,
            "phIncrement":0.5,
            "smiles": str(self.smiles),
            "username":"browser1",
            "elimAcid":[],
            "elimBase":[],
            "considerMethylAsAcid": True
        }
        post_string = json.dumps(sparc_post)
        logging.info("POST: {}".format(sparc_post))
        logging.info("POST (string): {}".format(post_string))
        try:
            response = requests.post(url, data=post_string, headers=headers, timeout=20)
            logging.info("response: {}".format(response.content))
            cleaned_json = response.content.replace(" ", "") # remove whitespace precaution..
            results = json.loads(cleaned_json)
        except Exception as e:
            logging.warning("SPARC PKA CALL ERROR: {}".format(e))
            raise
        else:
            return results


    def getPkaResults(self, results):
        """
        Gets pKa values from SPARC pKa request
        """
        try:
            pka_results = results['macroPkaResults'] # list of pka..
            pka_data = [] # return list of pka values..
            for pka in pka_results:
                if 'macroPka' in pka and pka['macroPka'] != -1000:
                    pka_data.append(pka['macroPka'])
                else:
                    pka_data.append("out of range")
        except Exception as e:
            logging.warning("Error getting pka from SPARC response: {}".format(e))
            raise Exception("error parsing sparc request")
        
        return {"pKa": pka_data}


    def makeCallForLogD(self):
        """
        Seprate call for octanol/water partition
        coefficient with pH (logD?)
        """
        logd_url = "/sparc-integration/rest/calc/logd"
        url = self.base_url + logd_url
        post = {
           "type":"LOGD",
           "solvent": {
              "solvents": None,
              "smiles": "OCCCCCCCC",
              "mixedSolvent": False,
              "name": "octanol"
           },
           "temperature": 25.0,
           "pH_minimum": 0,
           "pH_increment": 0.1,
           "ionic_strength": 0.0,
           "smiles": self.smiles
        }
        try:
            response = requests.post(url, data=json.dumps(post), headers=headers, timeout=20)
            results = json.loads(response.content)
            logging.info("{}".format(results))
        except Exception as e:
            logging.warning("SPARC LOGD CALL ERROR: {}".format(e))
            raise
        else:
            return results

    def getLogDForPH(self, results, ph=7.0):
        """
        Gets logD value at ph from
        logD response data
        TODO: add ph functionality
        """
        logging.info("getting sparc logd at ph: {}".format(ph))
        try:
            plot_data = results['plotCoordinates'] # list of [x,y]..
            logging.info("plot data: {}".format(plot_data))
            for xypair in plot_data:
                if xypair[0] == float(ph):
                    return xypair[1]
        except Exception as e:
            logging.warning("Error getting logD at PH from SPARC: {}".format(e))
            raise