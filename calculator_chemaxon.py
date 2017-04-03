import requests
import json
import logging
import os
import redis

from smilesfilter import parseSmilesByCalculator
from calculator import Calculator
import jchem_properties


# class JchemCalc(object):
class JchemCalc(Calculator):

    def __init__(self, prop_name=None):
        
        Calculator.__init__(self)  # inherit Calculator base class

        self.propsList = ['pKa', 'isoelectricPoint', 'majorMicrospecies', 'tautomerization', 'stereoisomer']  # speciation props
        
        self.baseUrl = os.environ['CTS_JCHEM_SERVER']
        self.name = ''
        self.max_retries = 3  # request retries if error
        self.timeout = 10  # request timeout in seconds
        self.structure = ''  # cas, smiles, iupac, etc.
        self.postData = {}
        self.results = ''
        self.methods = ['KLOP', 'VG', 'PHYS']  # kow_no_ph and kow_wph only
        self.props = ['water_sol', 'ion_con', 'kow_no_ph', 'kow_wph', 'water_sol_ph']  # available pchem props
        self.prop_name = prop_name  # prop name for JchemCalc instance

        # chemaxon speciation request
        self.speciation_request = {
            'run_type': None,
            'chem_struct': None,
            'smiles': None,
            'orig_smiles': None,
            'iupac': None,
            'formula': None,
            'mass': None,
            'get_pka': True,
            'get_taut': None,
            'get_stereo': None,
            'pKa_decimals': None,
            'pKa_pH_lower': None,
            'pKa_pH_upper': None,
            'pKa_pH_increment': None,
            'pH_microspecies': None,
            'isoelectricPoint_pH_increment': None,
            'tautomer_maxNoOfStructures': None,
            'tautomer_pH': None,
            'stereoisomers_maxNoOfStructures': None,
        }


    def setPostDataValue(self, propKey, propValue):
        """
        Can set one key:value with (propKey, propValue)
        """
        try:
            if propValue:
                self.postData[propKey] = propValue
        except KeyError as ke:
            logging.warning("key {} does not exist".format(propKey))
            raise
        except Exception as e:
            logging.warning("error occurred: {}".format(e))
            return None

    def setPostDataValues(self, multiKeyValueDict):
        """
        Can set multiple key:values at once w/ dict
        """
        try:
            for key, value in multiKeyValueDict.items():
                if value:
                    self.postData[key] = value
        except KeyError as ke:
            logging.warning("key {} does not exist".format(key))
            return None
        except Exception as e:
            logging.warning("error occurred: {}".format(e))
            return None

    def data_request_handler(self, request_dict, ws_protocol=False):
        """
        Handles requests to the JCHEM server.
          + For websocket protocol, it can handle a list of p-chem properties,
            and pushes the data up to the user.
          + For http protocol, it handles one p-chem property at a time and responds
            with a django.http HttpResponse.
        Inputs:
          + request_dict - POST data for p-chem request or speciation (transformation
            products moved to MetabolizerCalc)
          + ws_protocol - boolean for handling data for websocket or http
        """

        for key, val in self.pchem_request.items():
            if not key in request_dict.keys():
                logging.info("request key {} not in request, using default value: {}".format(key, val))
                request_dict.update({key: val})

        logging.info("Incoming request to ChemAxon's data_request_handler: {}".format(request_dict))

        _filtered_smiles = ''
        try:
            _filtered_smiles = parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
            # _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['smiles'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            request_dict.update({'data': 'Cannot filter SMILES for EPI data'})
            
            if not ws_protocol:
                return request_dict  # if not WS, just send object (for http/rest)

            self.redis_conn.publish(request_dict['sessionid'], json.loads(request_dict))
            return


        logging.info("Original CTS filtered SMILES: {}".format(request_dict['chemical']))
        logging.info("SMILES after filtering for ChemAxon calculator: {}".format(_filtered_smiles))


        if request_dict['service'] == 'getSpeciationData':

            data_obj = {
                'calc': "chemaxon", 
                'prop': "speciation_results",
                'node': request_dict['node'],
                'chemical': _filtered_smiles,
                'workflow': 'chemaxon',
                'run_type': 'batch'
            }

            # speciation data from chemspec model, for batch via ws
            # from cts_app.models.chemspec import chemspec_output

            # spec_inputs = request.POST.get('speciation_inputs')
            try:
                # todo: move error handling to modules that call calculator classes??? e.g., tasks.py, cts_rest.py
                _spec_inputs = request_dict['speciation_inputs']
            except KeyError as ke:
                logging.warning("speciation_inputs object needed for speciation request -- {}".format(ke))
                data_obj.update({'data': 'speciation POST needed'})
                # if not ws_protocol:
                #     return data_obj
                # self.redis_conn.publish(request_dict['sessionid'], json.dumps(data_obj))
                return data_obj

            _model_params = self.speciation_request
            for key, value in _model_params.items():
                if key in _spec_inputs:
                    _model_params.update({key: _spec_inputs[key]})

            _model_params.update({'smiles': _filtered_smiles, 'run_type': 'single'})

            logging.warning("MAKING REQUEST TO SPECIATION: {}".format(_model_params))

            # TODO: CTS API URL has env var somewhere...
            # chemspec_obj = chemspec_output.chemspecOutputPage(request)
            # speciation_response = requests.post('http://localhost:8000/cts/rest/speciation/run', data=json.dumps(_model_params))
            speciation_response = requests.post(
                                    os.environ.get('CTS_REST_SERVER') + '/cts/rest/speciation/run', 
                                    data=json.dumps(_model_params), 
                                    allow_redirects=True,
                                    verify=False)

            

            logging.warning("Speciation Response: {}".format(speciation_response.content))
            logging.warning("Response Type: {}".format(type(speciation_response.content)))

            speciation_data = json.loads(speciation_response.content)

            # logging.warning("SPECIATION DATA RECEIVED: {}".format(speciation_data))

            data_obj = {
                'calc': "chemaxon", 
                'prop': "speciation_results",
                'node': request_dict['node'],
                'chemical': _filtered_smiles,
                'workflow': 'chemaxon',
                'run_type': 'batch'
            }
            data_obj.update(speciation_data)

            # result_json = json.dumps(data_obj)
            # self.redis_conn.publish(request_dict['sessionid'], result_json)
            return data_obj

        else:

            # for prop in request_dict['props']:

            _response_dict = {}
            for key in request_dict.keys():
                _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

            _response_dict.update({'request_post': request_dict})
            # _response_dict.update({'prop': prop})

            logging.info("response dict: {}".format(_response_dict))

            try:

                if request_dict['prop'] == 'kow_wph' or request_dict['prop'] == 'kow_no_ph':
                    # for method in self.methods:

                    # request_dict.update({'method': method})
                    _response_dict.update({'method': method})
                    # _results = self.parse_jchem_results(request_dict)
                    # _results = self.getJchemPropData(request_dict['chemical'], prop, request_dict['ph'], method, request_dict['mass'])
                    _results = jchem_properties.JchemProperty().getJchemPropData(_response_dict)
                    logging.warning("RESULTS CALC CHEMAXON: {}".format(_results))

                    _response_dict.update({'data': _results['data'], 'method': method})

                    logging.warning("chemaxon _response_dict for {}: {}".format(request_dict['prop'], _response_dict))
                    # result_json = json.dumps(_response_dict)

                    # self.redis_conn.publish(request_dict['sessionid'], result_json)
                    return _response_dict


                else:
                    # _results = self.parse_jchem_results(request_dict)
                    # _results = self.getJchemPropData(request_dict['chemical'], prop, request_dict['ph'], None, request_dict['mass'])
                    _results = jchem_properties.JchemProperty().getJchemPropData(_response_dict)
                    _response_dict.update({'data': _results['data'], 'method': None})

                    logging.info("chemaxon _response_dict for {}: {}".format(request_dict['prop'], _response_dict))
                    return _response_dict
                    # result_json = json.dumps(_response_dict)
                    # self.redis_conn.publish(request_dict['sessionid'], result_json)

            except Exception as err:
                logging.warning("Exception occurred getting chemaxon data: {}".format(err))

                _response_dict.update({
                    'data': "cannot reach chemaxon calculator"
                })
                return _response_dict
                # result_json = json.dumps(_response_dict)
                # self.redis_conn.publish(request_dict['sessionid'], result_json)

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

    def booleanize(self, value):
        """
        django checkbox comes back as 'on' or 'off',
        or True/False depending on version, so this
        makes sure they're True/False
        """
        if value == 'on' or value == 'true':
            return True
        if value == 'off' or value == 'false':
            return False
        if isinstance(value, bool):
            return value

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

    # @classmethod
    # def getPropObject(self, prop):
    #     """
    #     For getting prop objects
    #     """
    #     if prop == 'pKa' or prop == 'ion_con':
    #         return jchem_properties.Pka()
    #     elif prop == 'isoelectricPoint':
    #         return jchem_properties.IsoelectricPoint()
    #     elif prop == 'majorMicrospecies':
    #         return jchem_properties.MajorMicrospecies()
    #     elif prop == 'tautomerization':
    #         return jchem_properties.Tautomerization()
    #     elif prop == 'stereoisomer':
    #         return jchem_properties.Stereoisomer()
    #     elif prop == 'solubility' or prop == 'water_sol' or prop == 'water_sol_ph':
    #         return jchem_properties.Solubility()
    #     elif prop == 'logP' or prop == 'kow_no_ph':
    #         return jchem_properties.LogP()
    #     elif prop == 'logD' or prop == 'kow_wph':
    #         return jchem_properties.LogD()
    #     else:
    #         raise ValueError("Error initializing jchem property class..")

    def getSpeciationResults(self, jchemResultObjects):
        """
        Loops jchemPropObjects (speciation results) from chemaxon,
        grabs the results and creates an object, jchemDictResults, that's
        used for chemspec_tables and data downloads.
        """
        jchem_results_obj = {}
        for key, value in jchemResultObjects.items():
            if value:
                if key == 'pKa':
                    jchem_results_obj.update({
                        'pka': jchemResultObjects['pKa'].getMostAcidicPka(),
                        'pkb': jchemResultObjects['pKa'].getMostBasicPka(),
                        'pka_parent': jchemResultObjects['pKa'].getParent(),
                        'pka_microspecies': jchemResultObjects['pKa'].getMicrospecies(),
                        'pka_chartdata': jchemResultObjects['pKa'].getChartData()
                    })
                elif key == 'majorMicrospecies':
                    jchem_results_obj.update({key: jchemResultObjects['majorMicrospecies'].getMajorMicrospecies()})
                elif key == 'isoelectricPoint':
                    jchem_results_obj.update({
                        key: jchemResultObjects['isoelectricPoint'].getIsoelectricPoint(),
                        'isopt_chartdata': jchemResultObjects['isoelectricPoint'].getChartData()
                    })
                elif key == 'tautomerization':
                    jchem_results_obj.update({'tautomers': jchemResultObjects['tautomerization'].getTautomers()})
                elif key == 'stereoisomers':
                    jchem_results_obj.update({key: jchemResultObjects['stereoisomers'].getStereoisomers()})

        # self.run_data.update(self.jchemDictResults)
        return jchem_results_obj
















# class JchemProperty(JchemCalc):
#     """
#     Jchem Physical-Chemical Properties
#     """
#     def __init__(self):


#         JchemCalc.__init__(self)
#         # self.request_timeout = 10
#         # self.headers = {'Content-Type': 'application/json'}
#         # self.max_retries = 3
#         # self.baseUrl = os.environ['CTS_JCHEM_SERVER']
#         self.url_pattern = '/webservices/rest-v0/util/calculate/{}'
#         self.results = ''  # json result
#         self.name = ''  # name of property
#         self.url = ''  # url to jchem ws endpoint
#         self.postData = {}  # POST data in json
#         self.ph = 7.0



#     def getJchemPropData(self, chemical, prop, ph=7.0, method=None, mass=None):
#         """
#         Calls jchem web services from chemaxon and
#         wraps data in a CTS data object (keys: calc, prop, method, data)
#         """

#         # resultDict = {"calc": "chemaxon", "prop": prop}

#         prop_obj = self.getPropObject(prop)
#         self.make_data_request(chemical, prop_obj, method)
#         prop_obj.get_data()

#         result_dict = {
#             'calc': 'chemaxon',
#             'prop': prop,
#             'data': result
#         }

#         # ADD METHOD KEY:VALUE IF LOGD OR LOGP...
#         if method:
#             result_dict['method'] = method

#         return result_dict


# class Pka(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'pKa'
#         self.url = self.url_pattern.format('pKa')
#         self.postData = {
#             "pHLower": 0.0,
#             "pHUpper": 14.0,
#             "pHStep": 0.1,
#             "temperature": 298.0,
#             "micro": False,
#             "considerTautomerization": True,
#             "pKaLowerLimit": 0.0,
#             "pKaUpperLimit": 14.0,
#             "prefix": "DYNAMIC"
#         }

#     def getMostAcidicPka(self):
#         """
#         Picks out pKa acidic value(s), returns list
#         """
#         pkaValList = []
#         if 'mostAcidic' in self.results:
#             for pkaVal in self.results['mostAcidic']:
#                 pkaValList.append(pkaVal)
#             return pkaValList
#         else:
#             logging.warning("key: 'mostAcidic' not in self.results")
#             return pkaValList

#     def getMostBasicPka(self):
#         """
#         Picks out pKa Basic value(s), returns list
#         """
#         pkaValList = []
#         if 'mostBasic' in self.results:
#             for pkaVal in self.results['mostBasic']:
#                 pkaValList.append(pkaVal)
#             return pkaValList
#         else:
#             logging.warning("no key 'mostBasic' in results")
#             return pkaValList

#     def getParent(self):
#         """
#         Gets parent image from result and adds structure
#         info such as formula, iupac, mass, and smiles.
#         Returns dict with keys: image, formula, iupac, mass, and smiles
#         """
#         try:
#             parentDict = {'image': self.results['result']['image']['image'], 'key': 'parent'}
#             parentDict.update(Calculator().getStructInfo(self.results['result']['structureData']['structure']))
#             return parentDict
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return None

#     def getMicrospecies(self):
#         """
#         Gets microspecies from pKa result and appends 
#         structure info (i.e., formula, iupac, mass, and smiles)
#         Returns list of microspeceies as dicts
#         with keys: image, formula, iupac, mass, and smiles
#         """
#         if 'microspecies' in self.results:
#             try:
#                 msList = []
#                 for ms in self.results['microspecies']:
#                     msStructDict = {}  # list element in msList
#                     msStructDict.update({'image': ms['image']['image'], 'key': ms['key']})
#                     structInfo = Calculator().getStructInfo(ms['structureData']['structure'])
#                     msStructDict.update(structInfo)
#                     msList.append(msStructDict)
#                 return msList
#             except KeyError as ke:
#                 logging.info("> key error: {}".format(ke))
#                 return None
#         else:
#             logging.info("no microspecies in results")
#             return None

#     def getChartData(self):
#         if 'chartData' in self.results:
#             microDistData = {}  # microspecies distribution data
#             for ms in self.results['chartData']:
#                 valuesList = []  # format: [[ph1,con1], [ph2, con2], ...] per ms
#                 for vals in ms['values']:
#                     xy = []  # [ph1, con1]
#                     xy.append(vals['pH'])
#                     xy.append(100.0 * vals['concentration'])  # convert to %
#                     valuesList.append(xy)
#                 microDistData.update({ms['key']: valuesList})
#             return microDistData
#         else:
#             return None


# class IsoelectricPoint(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'isoelectricPoint'
#         self.url = self.url_pattern.format('isoelectricPoint')
#         self.postData = {
#             "pHStep": 0.1,
#             "doublePrecision": 2
#         }

#     def getIsoelectricPoint(self):
#         """
#         Returns isoelectricPoint value from results
#         """
#         try:
#             return self.results['isoelectricPoint']
#         except KeyError:
#             logging.warning("key 'isoelectricPoint' not in results")
#             return None

#     def getChartData(self):
#         """
#         Returns isoelectricPoint chart data
#         """
#         valsList = []
#         try:
#             for pt in self.results['chartData']['values']:
#                 xyPair = []
#                 for key, value in pt.items():
#                     xyPair.append(value)
#                 valsList.append(xyPair)
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return valsList
#         else:
#             return valsList


# class MajorMicrospecies(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'majorMicrospecies'
#         self.url = self.url_pattern.format('majorMicrospecies')
#         self.postData = {
#             "pH": 7.0,
#             "takeMajorTautomericForm": False
#         }

#     def getMajorMicrospecies(self):
#         majorMsDict = {}
#         try:
#             majorMsDict.update({'image': self.results['result']['image']['image'], 'key': 'majorMS'})
#             structInfo = Calculator().getStructInfo(self.results['result']['structureData']['structure'])
#             majorMsDict.update(structInfo)  # add smiles, iupac, mass, formula key:values
#             return majorMsDict
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return None


# class Tautomerization(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'tautomerization'
#         self.url = self.url_pattern.format('tautomerization')
#         self.postData = {
#             "calculationType": "DOMINANT",
#             # "calculationType": "MAJOR",
#             "maxStructureCount": 1000,
#             "considerPH": False,
#             "enableMaxPathLength": True,
#             "maxPathLength": 4,
#             "rationalTautomerGenerationMode": False,
#             "singleFragmentMode": True,
#             "protectAromaticity": True,
#             "protectCharge": True,
#             "excludeAntiAromaticCompounds": True,
#             "protectDoubleBondStereo": False,
#             "protectAllTetrahedralStereoCenters": False,
#             "protectLabeledTetrahedralStereoCenters": False,
#             "protectEsterGroups": True,
#             "ringChainTautomerizationAllowed": False
#         }

#     def getTautomers(self):
#         """
#         returns dict w/ key 'tautStructs' and
#         value is a list of tautomer images
#         """
#         tautDict = {'tautStructs': [None]}
#         tautImageList = []
#         try:

#             tauts = self.results['result']  # for DOMINANT tautomers

#             for taut in tauts:
#                 tautStructDict = {'image': taut['image']['image'], 'key': 'taut'}
#                 structInfo = Calculator().getStructInfo(taut['structureData']['structure'])
#                 tautStructDict.update(structInfo)
#                 tautStructDict.update({'dist': 100 * round(taut['dominantTautomerDistribution'], 4)})
#                 tautImageList.append(tautStructDict)

#             tautDict.update({'tautStructs': tautImageList})
#             return tautImageList
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return None


# class Stereoisomer(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'stereoisomer'
#         self.url = self.url_pattern.format('stereoisomer')
#         self.postData = {
#             "stereoisomerismType": "TETRAHEDRAL",
#             "maxStructureCount": 100,
#             "protectDoubleBondStereo": False,
#             "protectTetrahedralStereo": False,
#             "filterInvalid3DStructures": False
#         }

#     def getStereoisomers(self):
#         stereoList = []
#         logging.warning("stereo results: {}".format(self.results))
#         try:
#             for stereo in self.results['result']:
#                 stereoDict = {'image': stereo['image']['image'], 'key': 'stereo'}
#                 structInfo = Calculator().getStructInfo(stereo['structureData']['structure'])
#                 stereoDict.update(structInfo)
#                 stereoList.append(stereoDict)
#             return stereoList
#         except KeyError as ke:
#             logging.warning("key error: {} @ jchem rest".format(ke))
#             return None


# class Solubility(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'solubility'
#         self.url = self.url_pattern.format('solubility')
#         self.postData = {
#             "pHLower": 0.0,
#             "pHUpper": 14.0,
#             "pHStep": 0.1,
#             "unit": "MGPERML"
#         }

#     def getIntrinsicSolubility(self):
#         """
#         Gets water solubility for chemaxon
#         """
#         try:
#             return 1000.0 * self.results['intrinsicSolubility']
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return None

#     def getPHDependentSolubility(self, ph=7.0):
#         """
#         Gets ph-dependent water solubility
#         """
#         try:
#             ws_list = self.results['pHDependentSolubility']['values']
#             ph = float(ph)
#             for item in ws_list:
#                 item_ph = item['pH']
#                 item_ws = item['solubility']
#                 if item_ph == ph:
#                     return item_ws
#             return "N/A"
#             # return 1000.0 * self.results['intrinsicSolubility']
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return None

#     def convertLogToMGPERL(self, log_val, mass):
#         """
#         Converts WS values of Log into mg/L.
#         Units of mass in mol/g
#         """
#         logging.info("MASS FOR CONVERSION: {}".format(mass))
#         logging.info("LOG VAL FOR CONVERSION: {}".format(log_val))
#         return 1000 * float(mass) * 10**(log_val)



# class LogP(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'logP'
#         self.url = self.url_pattern.format('logP')
#         self.methods = ['KLOP', 'VG', 'PHYS']
#         self.postData = {
#             "wVG": 1.0,
#             "wKLOP": 1.0,
#             "wPHYS": 1.0,
#             "Cl": 0.1,
#             "NaK": 0.1,
#             "considerTautomerization": False
#         }

#     def getLogP(self):
#         """
#         Gets pH-independent kow
#         """
#         try:
#             return self.results['logpnonionic']
#         except KeyError as ke:
#             logging.warning("ker error: {}".format(ke))
#             return None


# class LogD(JchemProperty):
#     def __init__(self):
#         JchemProperty.__init__(self)
#         self.name = 'logD'
#         self.url = self.url_pattern.format('logD')
#         self.methods = ['KLOP', 'VG', 'PHYS']
#         self.postData = {
#             "pHLower": 0.0,
#             "pHUpper": 14.0,
#             "pHStep": 0.1,
#             "wVG": 1.0,
#             "wKLOP": 1.0,
#             "wPHYS": 1.0,
#             "Cl": 0.1,
#             "NaK": 0.1,
#             "considerTautomerization": False
#         }

#     def getLogD(self, ph):
#         """
#         Gets pH-dependent kow
#         """
#         try:
#             ph = float(ph)
#             chartDataList = self.results['chartData']['values']
#             for xyPair in chartDataList:
#                 if xyPair['pH'] == round(ph, 1):
#                     value = xyPair['logD']
#                     break
#             return value
#         except KeyError as ke:
#             logging.warning("key error: {}".format(ke))
#             return None