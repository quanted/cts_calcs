__author__ = 'np'

"""
Classes for jchem ws - p-chem and other properties.
Intended to match the structure of the other
calculator classes.
"""

import requests
import json
import logging
import os
from jchem_rest import getStructInfo


# try:
#     from cts_app.cts_calcs.calculator import Calculator
# except ImportError as e:
#     from cts_calcs.calculator import Calculator


headers ={'Content-Type': 'application/json'}


class ChemaxonCalc(object):

    def __init__(self, prop_name=None):
        # speciation props:
        self.propsList = ['pKa', 'isoelectricPoint', 'majorMicrospecies', 'tautomerization', 'stereoisomer']
        
        self.baseUrl = os.environ['CTS_JCHEM_SERVER']
        self.name = ''
        self.url = ''
        self.max_retries = 3  # request retries if error
        self.timeout = 10  # request timeout in seconds
        self.structure = ''  # cas, smiles, iupac, etc.
        self.postData = {}
        self.results = ''
        self.propMap = {
            'water_sol': [],
            'ion_con': [], 
            'kow_no_ph': [],
            'kow_wph': []
        }
        self.props = ['water_sol', 'ion_con', 'kow_no_ph', 'kow_wph', 'water_sol_ph']  # available pchem props
        self.prop_name = prop_name  # prop name for ChemaxonCalc instance

        if prop_name in self.props:
            self.getPropObject(prop_name)  # __init__ automatically returns object
        #     return self.getPropObject(prop_name)

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

    def make_data_request(self, request_dict, message=None):
        """
        message - django channels message type. if message, send via ws/reply_channel attr
        """

        _response_dict = {}

         # populate p-chem response dict
        for key in request_dict.keys():
            _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

        _response_dict.update({'request_post': request_dict['request_post']})  # need this???


        for prop in request_dict['props']:

            try:

                if prop == 'kow_wph' or prop == 'kow_no_ph':
                    for method in self.methods:

                        request_dict.update({'method': method})
                        # _results = self.parse_jchem_results(request_dict)
                        _results = self.getJchemPropData(request_dict['chemical'], prop, request_dict['ph'], request_dict['mass'])

                        _response_dict.update({'data': _results['data'], 'method': method})

                        logging.info("chemaxon _results: {}".format(_results))
                        result_json = json.dumps(_response_dict)

                        self.redis_conn.publish(sessionid, result_json)
                        # if message:
                        #     message.reply_channel({'text': result_json})  # channels/ws
                        # else:
                        #     # can return because rest only handling one calc, prop, and method per request
                        #     return HttpResponse(result_json, content_type='application/json')

                else:
                    # _results = self.parse_jchem_results(request_dict)
                    _results = self.getJchemPropData(request_dict['chemical'], prop, request_dict['ph'], request_dict['mass'])
                    _response_dict.update({'data': _results['data']})

                    logging.info("chemaxon _results: {}".format(_results))

                    result_json = json.dumps(_response_dict)
                    self.redis_conn.publish(sessionid, result_json)
                    # if message:
                    #     message.reply_channel({'text': result_json})  # channels/ws
                    # else:
                    #     # can return because rest only handling one calc, prop, and method per request
                    #     return HttpResponse(result_json, content_type='application/json')

            except Exception as err:
                logging.warning("Exception occurred getting chemaxon data: {}".format(err))

                _response_dict.update({
                    'data': "cannot reach chemaxon calculator"
                })

                result_json = json.dumps(_response_dict)

                if message:
                    message.reply_channel({'text': result_json})  # channels/ws
                else:
                    # can return because rest only handling one calc, prop, and method per request
                    return HttpResponse(result_json, content_type='application/json')

    def makeDataRequest(self, structure, method=None, session=None):
        url = self.baseUrl + self.url
        self.postData.update({
            "result-display": {
                "include": ["structureData", "image"],
                "parameters": {
                    "structureData": "smiles"
                }
            }
        })
        postData = {
            "structure": structure,
            "parameters": self.postData
        }

        _valid_result = True  # for retry logic
        _retries = 0

        if method:
            postData['parameters']['method'] = method

        # retry data request to chemaxon server until max retries or a valid result is returned
        while not _valid_result or _retries < self.max_retries:
            try:
                response = requests.post(url, data=json.dumps(postData), headers=headers, timeout=self.timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    self.results = json.loads(response.content)
                    logging.info("Response from jchem server: {}".format(response))
                    break
                _retries += 1
            except Exception as e:
                logging.warning("Exception in jchem_calculator.py: {}".format(e))
                _retries += 1

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
            logging.warning("cts_celery jchem_calculator -- jchem server response status not 200, but: {}".format(response.status_code))
            # raise Exception("non 200 response from jchem: {}".format(response))
            return False
        
        # successful response, any further validating should go here (e.g., expected keys, error json from jchem server, etc.)
        # json_obj = json.loads(response.content)

        # TODO: verify if blank data, finding the source of the empty water sol values...
        return True

    @classmethod
    def getPropObject(self, prop):
        """
		For getting prop objects in a general,
		loop-friendly way
		"""
        if prop == 'pKa' or prop == 'ion_con':
            return Pka()
        elif prop == 'isoelectricPoint':
            return IsoelectricPoint()
        elif prop == 'majorMicrospecies':
            return MajorMicrospecies()
        elif prop == 'tautomerization':
            return Tautomerization()
        elif prop == 'stereoisomer':
            return Stereoisomer()
        elif prop == 'solubility' or prop == 'water_sol' or prop == 'water_sol_ph':
            return Solubility()
        elif prop == 'logP' or prop == 'kow_no_ph':
            return LogP()
        elif prop == 'logD' or prop == 'kow_wph':
            return LogD()
        else:
            raise ValueError("Error initializing jchem property class..")

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

    def getJchemPropData(self, chemical, prop, ph=7.0, method=None, sessionid=None, node=None, session=None, mass=None):
        """
        Calls jchem web services from chemaxon and
        wraps data in a CTS data object (keys: calc, prop, method, data)
        """

        resultDict = {"calc": "chemaxon", "prop": prop}

        result = ""
        if prop == 'water_sol':
            propObj = ChemaxonCalc.getPropObject('solubility')
            propObj.makeDataRequest(chemical, None, session)
            result = propObj.getIntrinsicSolubility()
        elif prop == 'ion_con':
            propObj = ChemaxonCalc.getPropObject('pKa')
            propObj.makeDataRequest(chemical, None, session)

            pkas = propObj.getMostAcidicPka() + propObj.getMostBasicPka()
            pkas.sort()

            result = {'pKa': pkas}
            # result = {'pKa': propObj.getMostAcidicPka(), 'pKb': propObj.getMostBasicPka()}
        elif prop == 'kow_no_ph':
            propObj = ChemaxonCalc.getPropObject('logP')
            propObj.makeDataRequest(chemical, method, session)
            result = propObj.getLogP()
        elif prop == 'kow_wph':
            propObj = ChemaxonCalc.getPropObject('logD')
            propObj.makeDataRequest(chemical, method, session)
            result = propObj.getLogD(ph)
        elif prop == 'water_sol_ph':
            propObj = ChemaxonCalc.getPropObject('solubility')
            propObj.makeDataRequest(chemical, method, session)
            result = propObj.getPHDependentSolubility(ph)
            result = propObj.convertLogToMGPERL(result, mass)
        else:
            result = None

        # ADD METHOD KEY:VALUE IF LOGD OR LOGP...
        resultDict['data'] = result
        if method:
            resultDict['method'] = method

        return resultDict

    def parse_jchem_results(self, request_dict):
        """
        Calls jchem web services from chemaxon and
        wraps data in a CTS data object (keys: calc, prop, method, data)
        """

        prop_obj = ChemaxonCalc.getPropObject(request_dict['prop'])
        url = self.baseUrl + prop_obj.url

        self.post_data.update({
            "result-display": {
                "include": ["structureData", "image"],
                "parameters": {
                    "structureData": "smiles"
                }
            }
        })
        _post_data = {
            'structure': request_dict['chemical'],
            'parameters': self.post_data
        }

        if 'method' in request_dict:
            _post_data['parameters']['method'] = request_dict['method']

        _valid_result = False  # for retry logic
        _retries = 0

        # retry data request to chemaxon server until max retries or a valid result is returned
        while not _valid_result or _retries < self.max_retries:
            try:
                response = requests.post(url, data=json.dumps(postData), headers=headers, timeout=self.timeout)
                _valid_result = self.validate_response(response)
                if _valid_result:
                    self.results = json.loads(response.content)
                    logging.info("Response from jchem server: {}".format(response))
                    break
                _retries += 1
            except Exception as e:
                logging.warning("Exception in jchem_calculator.py: {}".format(e))
                _retries += 1



class Pka(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'pKa'
        self.url = '/webservices/rest-v0/util/calculate/pKa'
        self.postData = {
            "pHLower": 0.0,
            "pHUpper": 14.0,
            "pHStep": 0.1,
            "temperature": 298.0,
            "micro": False,
            "considerTautomerization": True,
            "pKaLowerLimit": 0.0,
            "pKaUpperLimit": 14.0,
            "prefix": "DYNAMIC"
        }

    def getMostAcidicPka(self):
        """
		Picks out pKa acidic value(s), returns list
		"""
        pkaValList = []
        if 'mostAcidic' in self.results:
            # logging.info("$ type: {} $".format(self.results['mostAcidic']))
            for pkaVal in self.results['mostAcidic']:
                pkaValList.append(pkaVal)
            return pkaValList
        else:
            logging.warning("key: 'mostAcidic' not in self.results")
            return pkaValList

    def getMostBasicPka(self):
        """
		Picks out pKa Basic value(s), returns list
		"""
        pkaValList = []
        if 'mostBasic' in self.results:
            for pkaVal in self.results['mostBasic']:
                pkaValList.append(pkaVal)
            return pkaValList
        else:
            logging.warning("no key 'mostBasic' in results")
            return pkaValList

    def getParent(self):
        """
		Gets parent image from result and adds structure
		info such as formula, iupac, mass, and smiles.
		Returns dict with keys: image, formula, iupac, mass, and smiles
		"""
        try:
            parentDict = {'image': self.results['result']['image']['image'], 'key': 'parent'}
            parentDict.update(getStructInfo(self.results['result']['structureData']['structure']))
            return parentDict
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def getMicrospecies(self):
        """
		Gets microspecies from pKa result and appends 
		structure info (i.e., formula, iupac, mass, and smiles)
		Returns list of microspeceies as dicts
		with keys: image, formula, iupac, mass, and smiles
		"""
        if 'microspecies' in self.results:
            try:
                msList = []
                for ms in self.results['microspecies']:
                    msStructDict = {}  # list element in msList
                    msStructDict.update({'image': ms['image']['image'], 'key': ms['key']})
                    structInfo = getStructInfo(ms['structureData']['structure'])
                    msStructDict.update(structInfo)
                    msList.append(msStructDict)
                return msList
            except KeyError as ke:
                logging.info("> key error: {}".format(ke))
                return None
        else:
            logging.info("no microspecies in results")
            return None

    def getChartData(self):
        if 'chartData' in self.results:
            microDistData = {}  # microspecies distribution data
            for ms in self.results['chartData']:
                valuesList = []  # format: [[ph1,con1], [ph2, con2], ...] per ms
                for vals in ms['values']:
                    xy = []  # [ph1, con1]
                    xy.append(vals['pH'])
                    xy.append(100.0 * vals['concentration'])  # convert to %
                    valuesList.append(xy)
                microDistData.update({ms['key']: valuesList})
            return microDistData
        else:
            return None


class IsoelectricPoint(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'isoelectricPoint'
        self.url = '/webservices/rest-v0/util/calculate/isoelectricPoint'
        self.postData = {
            "pHStep": 0.1,
            "doublePrecision": 2
        }

    def getIsoelectricPoint(self):
        """
		Returns isoelectricPoint value from results
		"""
        try:
            return self.results['isoelectricPoint']
        except KeyError:
            logging.warning("key 'isoelectricPoint' not in results")
            return None

    def getChartData(self):
        """
		Returns isoelectricPoint chart data
		"""
        valsList = []
        try:
            for pt in self.results['chartData']['values']:
                xyPair = []
                for key, value in pt.items():
                    xyPair.append(value)
                valsList.append(xyPair)
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return valsList
        else:
            return valsList


class MajorMicrospecies(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'majorMicrospecies'
        self.url = '/webservices/rest-v0/util/calculate/majorMicrospecies'
        self.postData = {
            "pH": 7.0,
            "takeMajorTautomericForm": False
        }

    def getMajorMicrospecies(self):
        majorMsDict = {}
        try:
            majorMsDict.update({'image': self.results['result']['image']['image'], 'key': 'majorMS'})
            structInfo = getStructInfo(self.results['result']['structureData']['structure'])
            majorMsDict.update(structInfo)  # add smiles, iupac, mass, formula key:values
            return majorMsDict
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None


class Tautomerization(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'tautomerization'
        self.url = '/webservices/rest-v0/util/calculate/tautomerization'
        self.postData = {
            "calculationType": "DOMINANT",
            # "calculationType": "MAJOR",
            "maxStructureCount": 1000,
            "considerPH": False,
            "enableMaxPathLength": True,
            "maxPathLength": 4,
            "rationalTautomerGenerationMode": False,
            "singleFragmentMode": True,
            "protectAromaticity": True,
            "protectCharge": True,
            "excludeAntiAromaticCompounds": True,
            "protectDoubleBondStereo": False,
            "protectAllTetrahedralStereoCenters": False,
            "protectLabeledTetrahedralStereoCenters": False,
            "protectEsterGroups": True,
            "ringChainTautomerizationAllowed": False
        }

    def getTautomers(self):
        """
        returns dict w/ key 'tautStructs' and
        value is a list of tautomer images
        """
        tautDict = {'tautStructs': [None]}
        tautImageList = []
        try:
            # expecting list of result objects:
            # for taut in self.results['result']:
            # taut = self.results['result']  # single result object

            tauts = self.results['result']  # for DOMINANT tautomers

            # logging.warning("taut instance: {}".format(isinstance(taut, list)))

            # if isinstance(taut, list):
            #     taut = taut[0]
                
            # tautStructDict = {'image': taut['image']['image'], 'key': 'taut'}
            
            # structInfo = getStructInfo(taut['structureData']['structure'])
            # tautStructDict.update(structInfo)
            # # tautStructDict.update({'dist': 100 * round(taut['dominantTautomerDistribution'], 4)})
            # tautImageList.append(tautStructDict)


            # 10-19-16 request by Eric
            for taut in tauts:
                tautStructDict = {'image': taut['image']['image'], 'key': 'taut'}
                structInfo = getStructInfo(taut['structureData']['structure'])
                tautStructDict.update(structInfo)
                tautStructDict.update({'dist': 100 * round(taut['dominantTautomerDistribution'], 4)})
                tautImageList.append(tautStructDict)

            tautDict.update({'tautStructs': tautImageList})
            return tautImageList
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None


class Stereoisomer(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'stereoisomer'
        self.url = '/webservices/rest-v0/util/calculate/stereoisomer'
        self.postData = {
            "stereoisomerismType": "TETRAHEDRAL",
            "maxStructureCount": 100,
            "protectDoubleBondStereo": False,
            "protectTetrahedralStereo": False,
            "filterInvalid3DStructures": False
        }

    def getStereoisomers(self):
        stereoList = []
        try:
            for stereo in self.results['result']:
                stereoDict = {'image': stereo['image']['image'], 'key': 'stereo'}
                structInfo = getStructInfo(stereo['structureData']['structure'])
                stereoDict.update(structInfo)
                stereoList.append(stereoDict)
            return stereoList
        except KeyError as ke:
            logging.warning("key error: {} @ jchem rest".format(ke))
            return None


class Solubility(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'solubility'
        self.url = '/webservices/rest-v0/util/calculate/solubility'
        self.postData = {
            "pHLower": 0.0,
            "pHUpper": 14.0,
            "pHStep": 0.1,
            "unit": "MGPERML"
        }

    def getIntrinsicSolubility(self):
        """
		Gets water solubility for chemaxon
		"""
        try:
            logging.info("getting solubility from: {}".format(self.results))
            return 1000.0 * self.results['intrinsicSolubility']
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def getPHDependentSolubility(self, ph=7.0):
        """
        Gets ph-dependent water solubility
        """
        try:
            logging.info("getting solubility from: {}".format(self.results))
            ws_list = self.results['pHDependentSolubility']['values']
            ph = float(ph)
            for item in ws_list:
                item_ph = item['pH']
                item_ws = item['solubility']
                logging.info("ph: {}, ws: {}".format(item_ph, item_ws))
                # logging.info("types for ph: {}, ws: {}".format(type(item_ph), type(item_ws)))
                if item_ph == ph:
                    logging.info("getting solubility: {} at ph: {}".format(item_ws, item_ph))
                    return item_ws
            return "N/A"
            # return 1000.0 * self.results['intrinsicSolubility']
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None

    def convertLogToMGPERL(self, log_val, mass):
        """
        Converts WS values of Log into mg/L.
        Units of mass in mol/g
        """
        logging.info("MASS FOR CONVERSION: {}".format(mass))
        logging.info("LOG VAL FOR CONVERSION: {}".format(log_val))
        return 1000 * float(mass) * 10**(log_val)



class LogP(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'logP'
        self.url = '/webservices/rest-v0/util/calculate/logP'
        self.methods = ['KLOP', 'VG', 'PHYS']
        # logging.info("METHOD: {}".format(method))
        # if not method:
        #     logging.info("using WEIGHTED method for logP..")
        #     method = self.methods[3]
        # if not (method in self.methods): 
        #     key_err = "method {} not recognized as a method for logP ({}).."
        #     logging.warning(key_err.format(method, self.methods))
        #     raise KeyError(key_err.format(method, self.methods))
        self.postData = {
            "wVG": 1.0,
            "wKLOP": 1.0,
            "wPHYS": 1.0,
            "Cl": 0.1,
            "NaK": 0.1,
            "considerTautomerization": False
        }

    def getLogP(self):
        """
		Gets pH-independent kow
		"""
        try:
            return self.results['logpnonionic']
        except KeyError as ke:
            logging.warning("ker error: {}".format(ke))
            return None


class LogD(ChemaxonCalc):
    def __init__(self):
        ChemaxonCalc.__init__(self)
        self.name = 'logD'
        self.url = '/webservices/rest-v0/util/calculate/logD'
        self.methods = ['KLOP', 'VG', 'PHYS']
        self.postData = {
            "pHLower": 0.0,
            "pHUpper": 14.0,
            "pHStep": 0.1,
            "wVG": 1.0,
            "wKLOP": 1.0,
            "wPHYS": 1.0,
            "Cl": 0.1,
            "NaK": 0.1,
            "considerTautomerization": False
        }

    def getLogD(self, ph):
        """
		Gets pH-dependent kow
		"""
        try:
            ph = float(ph)
            chartDataList = self.results['chartData']['values']
            for xyPair in chartDataList:
                if xyPair['pH'] == round(ph, 1):
                    value = xyPair['logD']
                    break
            return value
        except KeyError as ke:
            logging.warning("key error: {}".format(ke))
            return None