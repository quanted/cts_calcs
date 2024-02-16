import requests
import json
import logging
import os
import redis
from .chemical_information import SMILESFilter
from .calculator import Calculator
from .jchem_properties import JchemProperty
from rdkit import Chem



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
        self.jchem_prop_obj = JchemProperty()

        # Chemaxon speciation request object:
        self.speciation_request = {
            'run_type': None,
            'chem_struct': None,
            'smiles': None,
            'orig_smiles': None,
            'iupac': None,
            'formula': None,
            'mass': None,
            'get_pka': None,
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


    def sort_microspecies(self, ms_list):
        """
        Sorts MS by pka using rdkit.
        """
        formalCharge={}
        for ms_obj in ms_list:
            chemical = ms_obj["structureData"]["structure"]
            smiles_resp = self.convertToSMILES({"chemical": chemical})
            smiles = smiles_resp.get("structure")
            mol = Chem.MolFromSmiles(smiles)
            fc = Chem.rdmolops.GetFormalCharge(mol)  # calculates formal charge
            ms_obj["fc"] = fc  # adding key:val for FC for each MS

        sorted_ms_list = sorted(ms_list, key=lambda item: item["fc"], reverse=True)
        sorted_ms_list = self.update_ms_id(sorted_ms_list)

        return sorted_ms_list

    def update_ms_id(self, ms_list):
        """
        """
        for i, ms in enumerate(ms_list):
            ms["key"] = "microspecies{}".format(i + 1)
        return ms_list 


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
                'run_type': "single",
                'request_post': request_dict
            }

            speciation_data = self.get_speciation_results(request_dict)

            data_obj['request_post'] = {'service': "speciation"}
            data_obj['data'] = speciation_data

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
                    _results = self.jchem_prop_obj.getJchemPropData(_response_dict)
                    _response_dict.update({'data': _results['data'], 'method': request_dict['method']})
                    return _response_dict

                else:
                    _results = self.jchem_prop_obj.getJchemPropData(_response_dict)
                    _response_dict.update({'data': _results['data'], 'method': None})
                    return _response_dict

            except Exception as err:
                logging.warning("Exception occurred getting chemaxon data: {}".format(err))

                _response_dict.update({
                    'data': "Cannot reach ChemAxon calculator"
                })
                return _response_dict



    def get_speciation_results(self, request):
        """
        Gets speciation results from jchem_properties.
        """
        jchemPropObjects = {}
        if 'speciation_inputs' in request:
            request.update(request['speciation_inputs'])
            del request['speciation_inputs']

        if request.get('get_pka'):
            # Makes call for pKa:
            pkaObj = JchemProperty.getPropObject('pKa')
            pkaObj.postData.update({
                "pHLower": request['pKa_pH_lower'],
                "pHUpper": request['pKa_pH_upper'],
                "pHStep": request['pKa_pH_increment'],
            })
            self.jchem_prop_obj.make_data_request(request['chemical'], pkaObj)
            jchemPropObjects['pKa'] = pkaObj

            ms = pkaObj.results["microspecies"]  # orig results, pre <img> wrappers and IDs
            sorted_ms_list = self.sort_microspecies(ms)  # sorts by FC
            sorted_ms_list = pkaObj.getMicrospecies()  # wraps MS for output page

            # Makes call for majorMS:
            majorMsObj = JchemProperty.getPropObject('majorMicrospecies')
            majorMsObj.postData.update({'pH': request['pH_microspecies']})
            self.jchem_prop_obj.make_data_request(request['chemical'], majorMsObj)
            jchemPropObjects['majorMicrospecies'] = majorMsObj

            # Makes call for isoPt:
            isoPtObj = JchemProperty.getPropObject('isoelectricPoint')
            isoPtObj.postData.update({'pHStep': request['isoelectricPoint_pH_increment']})
            self.jchem_prop_obj.make_data_request(request['chemical'], isoPtObj)
            jchemPropObjects['isoelectricPoint'] = isoPtObj

        if request.get('get_taut'):
            # Makes tautomer request:
            tautObj = JchemProperty.getPropObject('tautomerization')
            tautObj.postData.update({
                "maxStructureCount": request['tautomer_maxNoOfStructures'],
                "pH": request['tautomer_pH']
            })
            self.jchem_prop_obj.make_data_request(request['chemical'], tautObj)
            jchemPropObjects['tautomerization'] = tautObj

        if request.get('get_stereo'):
            # Makes stereoisomer request:
            stereoObj = JchemProperty.getPropObject('stereoisomer')
            stereoObj.postData.update({'maxStructureCount': request['stereoisomers_maxNoOfStructures']})
            self.jchem_prop_obj.make_data_request(request['smiles'], stereoObj)
            jchemPropObjects['stereoisomers'] = stereoObj

        speciation_results = self.jchem_prop_obj.getSpeciationResults(jchemPropObjects)

        return speciation_results