import json
import logging
import requests
import os
from calculator import Calculator
from calculator_measured import MeasuredCalc
from calculator_test import TestCalc
import smilesfilter


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
        self.props = ["water_sol", "vapor_press", "henrys_law_con", "mol_diss", "boiling_point"]
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
        self.request_timeout = 10

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

        # logging.info("Incoming request to SPARC's data_request_handler: {}".format(request_dict))

        _filtered_smiles = ''
        try:
            _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
            # _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['smiles'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            request_dict.update({'data': 'Cannot filter SMILES'})
            # self.redis_conn.publish(request_dict['sessionid'], json.loads(request_dict))
            return request_dict


        self.smiles = _filtered_smiles  # set smiles attribute to filtered smiles

        # Get melting point for sparc calculations.
        # Try Measured, then TEST..although it'll be slow
        melting_point = self.getMeltingPoint(_filtered_smiles, request_dict['sessionid'])
        # melting_point = 0.0  # TODO: add getMeltingPoint back after Measured and TEST refactor
        logging.warning("Using melting point: {} for SPARC calculation".format(melting_point))

        _response_dict = {}
        for key in request_dict.keys():
            _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

        _response_dict.update({'request_post': request_dict})
        # logging.info("response dict: {}".format(_response_dict))


        try:
            # if 'ion_con' in request_dict['props']:
            if request_dict.get('prop') == 'ion_con':
                response = self.makeCallForPka() # response as d ict returned..
                pka_data = self.getPkaResults(response)
                _response_dict.update({'data': pka_data, 'prop': 'ion_con'})
                return _response_dict

            # if 'kow_wph' in request_dict['props']:
            elif request_dict.get('prop') == 'kow_wph':
                response = self.makeCallForLogD() # response as dict returned..
                _response_dict.update({'data': self.getLogDForPH(response, request_dict['ph']), 'prop': 'kow_wph'})
                return _response_dict

            else:
                _post = self.get_sparc_query()
                _url = self.base_url + self.multiproperty_url
                _multi_response = self.makeDataRequest()

                if 'calculationResults' in _multi_response:
                    _multi_response = self.parseMultiPropResponse(_multi_response['calculationResults'], request_dict)
                    for prop_obj in _multi_response:
                        # if prop_obj['prop'] in request_dict['props'] and prop_obj['prop'] != 'ion_con' and prop_obj['prop'] != 'kow_wph':
                        if prop_obj['prop'] == request_dict['prop'] and prop_obj['prop'] != 'ion_con' and prop_obj['prop'] != 'kow_wph':
                            _prop = prop_obj['prop']
                            _data = prop_obj['data']
                            prop_obj.update(_response_dict)
                            prop_obj.update({'prop': _prop, 'data': _data})
                            return prop_obj

        except Exception as err:
            logging.warning("!!! Exception occurred getting SPARC data: {} !!!".format(err))

            _response_dict.update({
                'data': "request timed out",
                'prop': request_dict.get('prop')
            })

            return _response_dict


    def makeDataRequest(self):
        _post = self.get_sparc_query()
        _url = self.base_url + self.multiproperty_url

        logging.info("SPARC URL: {}".format(_url))
        logging.info("SPARC POST: {}".format(_post))

        return self.request_logic(_url, _post)


    def request_logic(self, url, post_data):
        """
        Handles retries and validation of responses
        """

        _valid_result = False  # for retry logic
        _retries = 0
        while not _valid_result and _retries < self.max_retries:
            # retry data request to chemaxon server until max retries or a valid result is returned
            try:
                response = requests.post(url, data=json.dumps(post_data), headers=self.headers, timeout=self.request_timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    self.results = json.loads(response.content)
                    # break
                    return self.results
                _retries += 1
            except Exception as e:
                logging.warning("Exception in calculator_sparc.py: {}".format(e))
                _retries += 1

            logging.info("Max retries: {}, Retries left: {}".format(self.max_retries, _retries))
        self.results = "calc server not found"
        return self.results


    def validate_response(self, response):
        """
        Validates sparc response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("sparc server response status: {}".format(response.status_code))
            return False
        
        try:
            response_obj = json.loads(response.content)
        except Exception as e:
            logging.warning("Could not convert response to json object, sparc validate_response: {}".format(e))
            return False

        response_type = response_obj.get('type')  # get SPARC property name

        logging.info("Validating {} property".format(response_type))

        # prop specific validation:
        if response_type == 'LOGD':
            # check 'plotCoordinates' key, should be list and not None
            if not isinstance(response_obj.get('plotCoordinates'), list):
                # retry if no 'plotCoordinates' key or 'plotCoordinates' is None
                logging.warning("SPARC LOGD 'plotCoordinates' not list as expected...Retrying request...")
                return False

        # successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
        # json_obj = json.loads(response.content)

        # TODO: verify if blank data, finding the source of the empty water sol values...
        return True


    def parseMultiPropResponse(self, results, request_dict):
        """
        Loops through data grabbing the results
        and building {calc, prop, data} objects
        for front end.

        TODO: Add more info to returned data (object instead of value)
        """
        if not results or not isinstance(results, list):
            raise Exception("(sparc) no results given to parse or not a list..")

        sparc_response = []
        sparc_response_props = []  # list of multi response props to make sure all props made it
        # logging.info("parsing results: {}".format(results))
        for prop_data in results:
            # logging.info("PROP DATA: {}".format(prop_data))
            sparc_prop = prop_data['type']
            logging.info("sparc prop: {}".format(sparc_prop))
            cts_prop_name = self.sparc_props.get(sparc_prop)
            data = prop_data['result']
            logging.info("cts prop name: {}".format(cts_prop_name))
            if cts_prop_name:
                data_obj = {'calc': 'sparc', 'prop': cts_prop_name, 'data': data}
                logging.info("data obj: {}".format(data_obj))
                sparc_response.append(data_obj)
                sparc_response_props.append(cts_prop_name)
                # logging.info("sparc response: {}".format(sparc_response))

        # logging.info("sparc multi response props list: {}".format(sparc_response_props))
        # logging.info("user request props list: {}".format(request_dict['props']))

        for prop in request_dict['props']:
            if not prop in sparc_response_props and prop != 'ion_con' and prop != 'kow_wph':
                # if sparc response doesn't have user request prop from multi-response, PANIC!
                logging.info("requested prop {} missing from sparc multi response...".format(prop))

                # add data obj for missing prop with error message for user:
                data_obj = {'calc': 'sparc', 'prop': prop, 'data': "prop not found"}
                sparc_response.append(data_obj)

        # logging.info("SPARC RESPONSE: {}".format(sparc_response))
        return sparc_response


    def makeCallForPka(self):
        """
        Separate call for SPARC pKa
        """
        _pka_url = "/sparc-integration/rest/calc/fullSpeciation"
        _url = self.base_url + _pka_url
        logging.info("URL: {}".format(_url))
        _sparc_post = {
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
        _post_string = json.dumps(_sparc_post)

        return self.request_logic(_url, _sparc_post)


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
        _logd_url = "/sparc-integration/rest/calc/logd"
        _url = self.base_url + _logd_url
        _post = {
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

        logd_results = self.request_logic(_url, _post)
        # logging.info("LOGD RESULTs: {}".format(logd_results))
        return logd_results


    def getLogDForPH(self, results, ph=7.0):
        """
        Gets logD value at ph from
        logD response data
        TODO: add ph functionality
        """
        # logging.info("getting sparc logd at ph: {}".format(ph))
        try:
            plot_data = results['plotCoordinates'] # list of [x,y]..
            # logging.info("plot data: {}".format(plot_data))
            for xypair in plot_data:
                if xypair[0] == float(ph):
                    return xypair[1]
        except Exception as e:
            logging.warning("LOGD ERROR RESULTS: {}".format(results))
            logging.warning("Error getting logD at PH from SPARC: {}".format(e))
            raise

    def getMeltingPoint(self, structure, sessionid):
        """
        Gets mass of structure from Measured, tries
        TEST if not available in Measured. Returns 0.0
        if neither have mp value.
        """
        melting_point_request = {
            'calc': "measured",  # should prob be measured
            # 'props': ['melting_point'],
            'prop': 'melting_point',
            'chemical': structure,
            'sessionid': sessionid
        }
        # todo: catch measured errors, then try epi melting point..
        # request = NotDjangoRequest(melting_point_request)
        # melting_point_response = measured_views.request_manager(request)
        measured_mp_response = MeasuredCalc().data_request_handler(melting_point_request)

        # # convert to python dict
        try:
            melting_point = json.loads(measured_mp_response.content)['data']
        except Exception as e:
            logging.warning("Error in sparc_cts/worker.py: {}".format(e))
            melting_point = 0.0

        logging.warning("MELTING POINT RESPONSE: {}".format(measured_mp_response))
        logging.warning("MELTING POINT RESPONSE TYPE: {}".format(type(measured_mp_response)))

        if not isinstance(melting_point, float):
            logging.warning("Trying to get MP from TEST..")
            try:
                melting_point_request['calc'] = 'test'
                # request = NotDjangoRequest(melting_point_request)
                # test_melting_point_response = test_views.request_manager(request)
                test_mp_response = TestCalc().data_request_handler(melting_point_request)
                logging.warning("TEST MP RESPONSE CONTENT: {}".format(test_melting_point_response.content))
                melting_point = json.loads(test_melting_point_response.content)[0]['data']
                logging.warning("TEST MP VALUE: {}".format(melting_point))
            except Exception as e:
                logging.warning("Error in sparc_cts/worker.py: {}".format(e))
                melting_point = 0.0

            logging.warning("TEST MP TYPE: {}:".format(type(melting_point)))

            if not isinstance(melting_point, float):
                melting_point = 0.0
        # else:
        #     melting_point = melting_point_obj['data']

        logging.warning("MELTING POINT VALUE: {}".format(melting_point))

        return melting_point