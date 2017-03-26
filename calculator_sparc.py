import json
import logging
import requests
import os
from calculator import Calculator


class SparcCalc(Calculator):
    def __init__(self, smiles=None, meltingpoint=0.0, pressure=760.0, temperature=25.0):

        Calculator.__init__(self)  # inherit Calculator base class

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


    def data_request_handler(self, request_dict, message=None):

        for key, val in self.pchem_request.items():
            if not key in request_dict.keys():
                logging.info("request key {} not in request, using default value: {}".format(key, val))
                request_dict.update({key: val})

        logging.info("Incoming request to ChemAxon's data_request_handler: {}".format(request_dict))

        _filtered_smiles = ''
        try:
            _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            request_dict.update({'data': 'Cannot filter SMILES for EPI data'})
            self.redis_conn.publish(request_dict['sessionid'], json.loads(request_dict))


        # Get melting point for sparc calculations.
        # Try Measured, then TEST..although it'll be slow
        melting_point = getMass(_filtered_smiles, request_dict['sessionid'])
        logging.warning("Using melting point: {} for SPARC calculation".format(melting_point))

        calcObj = SparcCalc(_filtered_smiles, meltingpoint=melting_point)
        
        multi_response = None
        sparc_results = []

        # synchronous calls for ion_con, kow_wph, and the rest:
        # don't need this loop, just do "if 'ion_con' in prop: make request"

        logging.warning("sparc props: {}".format(props))

        # sparc_response = {
        #     'calc': 'sparc',
        #     # 'prop': prop,
        #     'node': node,
        #     'chemical': structure,
        #     'request_post': request.POST,
        #     'run_type': run_type,
        #     'workflow': workflow
        # }

        if run_type == 'rest':
            props = [prop]



        _response_dict = {}
        for key in request_dict.keys():
            _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

        _response_dict.update({'request_post': request_dict})
        _response_dict.update({'prop': prop})

        logging.info("response dict: {}".format(_response_dict))


        try:
            if 'ion_con' in props:
                response = calcObj.makeCallForPka() # response as d ict returned..
                pka_data = calcObj.getPkaResults(response)
                _response_dict.update({'data': pka_data, 'prop': 'ion_con'})
                logging.info("ion_con response: {}".format(_response_dict))
                result_json = json.dumps(_response_dict)
                redis_conn.publish(request_dict['sessionid'], result_json)

            if 'kow_wph' in props:
                ph = request.POST.get('ph') # get PH for logD calculation..
                response = calcObj.makeCallForLogD() # response as dict returned..
                _response_dict.update({'data': calcObj.getLogDForPH(response, ph), 'prop': 'kow_wph'})
                logging.info("kow_wph response: {}".format(_response_dict))
                result_json = json.dumps(_response_dict)
                redis_conn.publish(request_dict['sessionid'], result_json)

            multi_response = calcObj.makeDataRequest()
            logging.info("MULTI RESPONSE: {}".format(multi_response))
            if 'calculationResults' in multi_response:
                multi_response = calcObj.parseMultiPropResponse(multi_response['calculationResults'])
                logging.info("Parsed Multi Response: {}".format(multi_response))
                for prop_obj in multi_response:
                    logging.info("PROP OBJECT: {}".format(prop_obj))
                    logging.info("requested props: {}".format(props))
                    logging.info("prop in props: {}".format(prop_obj['prop'] in props))
                    logging.info("prop: {}".format(prop_obj['prop']))
                    if prop_obj['prop'] in props and prop_obj['prop'] != 'ion_con' and prop_obj['prop'] != 'kow_wph':
                        prop_obj.update({'node': node, 'chemical': structure, 'request_post': request.POST, 'workflow': workflow, 'run_type': run_type})
                        # prop_obj.update(sparc_response)
                        logging.info("multiprop response: {}".format(prop_obj))
                        result_json = json.dumps(prop_obj) 
                        redis_conn.publish(request_dict['sessionid'], result_json)

        except Exception as err:
            logging.warning("Exception occurred getting SPARC data: {}".format(err))

            for prop in props:

                post_data.update({
                    'data': "data request timed out",
                    'prop': prop,
                    'request_post': request.POST
                })

                redis_conn.publish(request_dict['sessionid'], json.dumps(post_data))


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