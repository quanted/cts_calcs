import requests
import json
import logging
import os
import redis

import smilesfilter
import jchem_properties
import jchem_rest
import data_walks
from calculator import Calculator


# class ChemaxonCalc(object):
class ChemaxonCalc(Calculator):

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
        self.prop_name = prop_name  # prop name for ChemaxonCalc instance

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

        self.transformation_request = {
            'structure': None,
            'generationLimit': 1,  # make sure to get this from front end
            'populationLimit': 0,
            'likelyLimit': 0.001,
            'transformationLibraries': ["hydrolysis", "abiotic_reduction"],  # NOTE: no transformationLibraries key:val for mammalian metabolism
            'excludeCondition': ""
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

    def data_request_handler(self, request_dict, message=None):
        """
        message - django channels message type. if message, send via ws/reply_channel attr
        """

        for key, val in self.pchem_request.items():
            if not key in request_dict.keys():
                logging.info("request key {} not in request, using default value: {}".format(key, val))
                request_dict.update({key: val})

        logging.info("Incoming request to ChemAxon's data_request_handler: {}".format(request_dict))

        _filtered_smiles = ''
        try:
            _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
            # _filtered_smiles = smilesfilter.parseSmilesByCalculator(request_dict['smiles'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            request_dict.update({'data': 'Cannot filter SMILES for EPI data'})
            self.redis_conn.publish(request_dict['sessionid'], json.loads(request_dict))


        logging.info("Original CTS filtered SMILES: {}".format(request_dict['chemical']))
        logging.info("SMILES after filtering for ChemAxon calculator: {}".format(_filtered_smiles))

        # request_dict['chemical'] = _filtered_smiles  # use calc-filtered smiles for cts-filtered for all calcs in output????

        if request_dict['service'] == 'getTransProducts':
            # getTransProducts chemaxon service via ws..
            request = {
                'structure': _filtered_smiles,
                'generationLimit': 1,  # make sure to get this from front end
                'populationLimit': 0,
                'likelyLimit': 0.001,
                # 'transformationLibraries': ['human_biotransformation'],  # TODO: Use UI lib choices
                'excludeCondition': ""  # 'generateImages': False
            }

            response = jchem_rest.getTransProducts(request)
            data_walks.j = 0
            data_walks.metID = 0
            results = data_walks.recursive(response, 1)

            data_obj = {
                'calc': "chemaxon", 
                'prop': "products",
                'node': request_dict['node'],
                'data': json.loads(results),
                'chemical': _filtered_smiles,
                'workflow': 'gentrans',
                'run_type': 'batch'
            }
            
            result_json = json.dumps(data_obj)
            logging.info("publishing to redis")
            self.redis_conn.publish(request_dict['sessionid'], result_json)

        elif request_dict['service'] == 'getSpeciationData':

            # TODO: Get this through jchem_rest and not the chemspec model!!!

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
                self.redis_conn.publish(request_dict['sessionid'], json.dumps(data_obj))

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

            result_json = json.dumps(data_obj)
            self.redis_conn.publish(request_dict['sessionid'], result_json)

        else:

            for prop in request_dict['props']:

                _response_dict = {}
                for key in request_dict.keys():
                    _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

                _response_dict.update({'request_post': request_dict})
                _response_dict.update({'prop': prop})

                logging.info("response dict: {}".format(_response_dict))

                try:

                    if prop == 'kow_wph' or prop == 'kow_no_ph':
                        for method in self.methods:

                            request_dict.update({'method': method})
                            # _results = self.parse_jchem_results(request_dict)
                            _results = self.getJchemPropData(request_dict['chemical'], prop, request_dict['ph'], method, request_dict['mass'])

                            _response_dict.update({'data': _results['data'], 'method': method})

                            logging.info("chemaxon _response_dict for {}: {}".format(prop, _response_dict))
                            result_json = json.dumps(_response_dict)

                            self.redis_conn.publish(request_dict['sessionid'], result_json)
                            # if message:
                            #     message.reply_channel({'text': result_json})  # channels/ws
                            # else:
                            #     # can return because rest only handling one calc, prop, and method per request
                            #     return HttpResponse(result_json, content_type='application/json')

                    else:
                        # _results = self.parse_jchem_results(request_dict)
                        _results = self.getJchemPropData(request_dict['chemical'], prop, request_dict['ph'], None, request_dict['mass'])
                        _response_dict.update({'data': _results['data'], 'method': None})

                        logging.info("chemaxon _response_dict for {}: {}".format(prop, _response_dict))

                        result_json = json.dumps(_response_dict)
                        self.redis_conn.publish(request_dict['sessionid'], result_json)
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
                    self.redis_conn.publish(request_dict['sessionid'], result_json)
                    # if message:
                    #     message.reply_channel({'text': result_json})  # channels/ws
                    # else:
                    #     # can return because rest only handling one calc, prop, and method per request
                    #     return HttpResponse(result_json, content_type='application/json')

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

    @classmethod
    def getPropObject(self, prop):
        """
        For getting prop objects
        """
        if prop == 'pKa' or prop == 'ion_con':
            return jchem_properties.Pka()
        elif prop == 'isoelectricPoint':
            return jchem_properties.IsoelectricPoint()
        elif prop == 'majorMicrospecies':
            return jchem_properties.MajorMicrospecies()
        elif prop == 'tautomerization':
            return jchem_properties.Tautomerization()
        elif prop == 'stereoisomer':
            return jchem_properties.Stereoisomer()
        elif prop == 'solubility' or prop == 'water_sol' or prop == 'water_sol_ph':
            return jchem_properties.Solubility()
        elif prop == 'logP' or prop == 'kow_no_ph':
            return jchem_properties.LogP()
        elif prop == 'logD' or prop == 'kow_wph':
            return jchem_properties.LogD()
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

    def getJchemPropData(self, chemical, prop, ph=7.0, method=None, mass=None):
        """
        Calls jchem web services from chemaxon and
        wraps data in a CTS data object (keys: calc, prop, method, data)
        """

        resultDict = {"calc": "chemaxon", "prop": prop}

        result = ""
        if prop == 'water_sol':
            propObj = ChemaxonCalc.getPropObject('solubility')
            self.make_data_request(chemical, propObj, None)
            result = propObj.getIntrinsicSolubility()
        elif prop == 'ion_con':
            propObj = ChemaxonCalc.getPropObject('pKa')
            self.make_data_request(chemical, propObj, None)

            pkas = propObj.getMostAcidicPka() + propObj.getMostBasicPka()
            pkas.sort()

            result = {'pKa': pkas}
            # result = {'pKa': propObj.getMostAcidicPka(), 'pKb': propObj.getMostBasicPka()}
        elif prop == 'kow_no_ph':
            propObj = ChemaxonCalc.getPropObject('logP')
            self.make_data_request(chemical, propObj, method)
            result = propObj.getLogP()
        elif prop == 'kow_wph':
            propObj = ChemaxonCalc.getPropObject('logD')
            self.make_data_request(chemical, propObj, method)
            result = propObj.getLogD(ph)
        elif prop == 'water_sol_ph':
            propObj = ChemaxonCalc.getPropObject('solubility')
            self.make_data_request(chemical, propObj, None)
            result = propObj.getPHDependentSolubility(ph)
            result = propObj.convertLogToMGPERL(result, mass)
        else:
            result = None

        # ADD METHOD KEY:VALUE IF LOGD OR LOGP...
        resultDict['data'] = result
        if method:
            resultDict['method'] = method

        return resultDict