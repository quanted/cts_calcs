import requests
import json
import logging
import os
import redis
from .chemical_information import SMILESFilter
from .calculator import Calculator
from .jchem_properties import JchemProperty



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
        self.format_url = '/rest-v0/util/analyze'  # returns chemical's format (e.g., "smiles", "casrn")

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



    def data_request_handler(self, request_dict):
        """
        Handles requests to the JCHEM server.
        Inputs:
          + request_dict - POST data for p-chem request or speciation (transformation
            products moved to MetabolizerCalc).
        """

        for key, val in self.pchem_request.items():
            if not key in request_dict.keys():
                request_dict.update({key: val})

        _filtered_smiles = ''
        try:
            _filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            request_dict.update({'data': 'Cannot filter SMILES for ChemAxon data'})
            return request_dict  # if not WS, just send object (for http/rest)

        if request_dict['service'] == 'getSpeciationData':

            data_obj = {
                'calc': "chemaxon", 
                'prop': "speciation_results",
                'node': request_dict['node'],
                'chemical': _filtered_smiles,
                'workflow': 'chemaxon',
                'run_type': 'batch'
            }

            try:
                _spec_inputs = request_dict['speciation_inputs']
            except KeyError as ke:
                logging.warning("speciation_inputs object needed for speciation request -- {}".format(ke))
                data_obj.update({'data': 'speciation POST needed'})
                return data_obj

            _model_params = self.speciation_request
            for key, value in _model_params.items():
                if key in _spec_inputs:
                    _model_params.update({key: _spec_inputs[key]})

            _model_params.update({
                'smiles': _filtered_smiles,
                'run_type': 'single',
                'chemical': request_dict['chemical'],
                'method': 'POST',
            })


            # TODO: CTS API URL has env var somewhere...
            # chemspec_obj = chemspec_output.chemspecOutputPage(request)
            # speciation_response = requests.post('http://localhost:8000/cts/rest/speciation/run', data=json.dumps(_model_params))
            speciation_response = requests.post(
                                    os.environ.get('CTS_REST_SERVER') + '/cts/rest/speciation/run', 
                                    data=json.dumps(_model_params), 
                                    allow_redirects=True,
                                    verify=False)


            speciation_data = json.loads(speciation_response.content)

            data_obj = {
                'calc': "chemaxon", 
                'prop': "speciation_results",
                'node': request_dict['node'],
                'chemical': _filtered_smiles,
                'workflow': 'chemaxon',
                'run_type': 'batch',
                'request_post': {'service': "speciation"}
            }
            data_obj.update(speciation_data)

            return data_obj

        else:

            _response_dict = {}
            for key in request_dict.keys():
                if not key == 'nodes':
                    _response_dict[key] = request_dict.get(key)  # fill any overlapping keys from request1

            _response_dict.update({'request_post': request_dict})

            try:

                if request_dict['prop'] == 'kow_wph' or request_dict['prop'] == 'kow_no_ph':
                    _response_dict.update({'method': request_dict['method']})
                    _results = JchemProperty().getJchemPropData(_response_dict)
                    _response_dict.update({'data': _results['data'], 'method': request_dict['method']})
                    return _response_dict

                else:
                    _results = JchemProperty().getJchemPropData(_response_dict)
                    _response_dict.update({'data': _results['data'], 'method': None})
                    return _response_dict

            except Exception as err:
                logging.warning("Exception occurred getting chemaxon data: {}".format(err))

                _response_dict.update({
                    'data': "Cannot reach ChemAxon calculator"
                })
                return _response_dict



    def validate_response(self, response):
        """
        Validates jchem response.
        Returns False if data is null, or any other
        values that indicate an error
        """
        if response.status_code != 200:
            logging.warning("cts_celery calculator_chemaxon -- jchem server response status not 200, but: {}".format(response.status_code))
            return False

        # TODO: verify if blank data, finding the source of the empty water sol values... 
        return True