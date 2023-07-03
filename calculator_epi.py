import requests
import json
import logging
import os
from rdkit import Chem

from .calculator import Calculator
from .chemical_information import SMILESFilter



class RdkitCalc:
    """
    rdkit
    """
    def __init__(self):
        # Establish smarts objects, NOTE these strings will not change!
        self.CAE_smarts = '[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]'
        self.CAE = Chem.MolFromSmarts(self.CAE_smarts)
        self.CarbAnhydride_smarts = '[CX3;$([H0][#6]),$([H1])](=[OX1])[#8X2][CX3;$([H0][#6]),$([H1])](=[OX1])'
        self.CarbAnhydride = Chem.MolFromSmarts(self.CarbAnhydride_smarts)

    def get_functional_groups_anhydride(self, smiles):
        """
        List of atoms in anhydride functional group.
        """
        mol = Chem.MolFromSmiles(smiles)
        anhydride_atom = list(Chem.Mol.GetSubstructMatches(mol, self.CarbAnhydride, uniquify=True)[0])
        logging.info("Anhydride group: {}".format(anhydride_atom))
        return anhydride_atom

    def get_functional_groups_cae(self, smiles):
        """
        List of atoms in carboxylic acid ester functional group.
        """
        mol = Chem.MolFromSmiles(smiles)
        CAE_atom = list(Chem.Mol.GetSubstructMatches(mol, self.CAE, uniquify=True)[0])
        logging.info("CAE group: {}".format(CAE_atom))
        return CAE_atom

    def get_functional_groups(self, route, smiles):
        """
        Calls rdkit to get functional groups.
        """
        func_group = []
        if "anhydride" in route.lower():
            func_group = self.get_functional_groups_anhydride(smiles)
        elif route.lower() == "carboxylic acid ester hydrolysis":
            func_group = self.get_functional_groups_cae(smiles)
        return func_group


class EpiCalc(Calculator):
    """
    EPI Suite Calculator
    """
    def __init__(self):
        Calculator.__init__(self)
        self.rdkit = RdkitCalc()
        self.method = None
        self.postData = {"smiles" : ""}
        self.name = "epi"
        self.baseUrl = os.environ['CTS_EPI_SERVER']
        # self.urlStruct = "/episuiteapi/rest/episuite/estimated"  # newest way - server
        # self.urlStruct = "/rest/episuite/estimated"  # newest way - local
        self.methods = None
        self.melting_point = None
        self.epi_props = ['melting_point', 'boiling_point', 'water_solubility', 'vapor_pressure', 'henrys_law_constant', 'log_kow', 'log_koc', 'log_bcf', 'log_baf']
        self.props = ['melting_point', 'boiling_point', 'water_sol', 'vapor_press', 'henrys_law_con', 'kow_no_ph', 'koc', 'log_bcf', 'log_baf']
        self.propMap = {
            'melting_point': {
               'result_key': 'melting_point'
            },
            'boiling_point': {
               'result_key': 'boiling_point'
            },
            'water_sol': {
               'result_key': 'water_solubility',
               'methods': {'WSKOW': "WSKOW", 'WATERNT': "WATERNT"}
            },
            'vapor_press': {
               'result_key': 'vapor_pressure'
            },
            'henrys_law_con': {
                'result_key': 'henrys_law_constant'
            },
            'kow_no_ph': {
                'result_key': 'log_kow'
            },
            'koc': {
                'result_key': 'log_koc',
                'methods': {'MCI': "MCI", 'Kow': "KOW"}
            },
            'log_bcf': {
                'result_key': 'log_bcf',
                'methods': {'regression': "REG", 'Arnot-Gobas': "A-G"}
            },
            'log_baf': {
                'result_key': 'log_baf',
                'methods': {'Arnot-Gobas': "A-G"}
            },
            'qsar': {
                'result_key': 'qsar',
            }
        }
        self.qsar_request_map = {
            'halogenated aliphatics: elimination': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (no adjacent x)': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (vicinal x)': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (geminal x)': 'hydrolysis/alkylhalide',
            'epoxide hydrolysis': 'hydrolysis/epoxide',
            'organophosphorus ester hydrolysis 1': 'hydrolysis/phosphate',  # Phosphate or Thiophosphate
            'organophosphorus ester hydrolysis 2': 'hydrolysis/phosphate',  # Phosphate or Thiophosphate
            'carboxylic acid ester hydrolysis': 'hydrolysis/ester',
            'anhydride hydrolysis': 'hydrolysis/anhydride',
            'carbamate hydrolysis': 'hydrolysis/carbamate'
        }
        self.cleaved_list = [
            'organophosphorus ester hydrolysis 1',
            'organophosphorus ester hydrolysis 2',
            'carboxylic acid ester hydrolysis',
            'anhydride hydrolysis',
            'carbamate hydrolysis'
        ]
        self.op_esters = [
            'organophosphorus ester hydrolysis 1',
            'organophosphorus ester hydrolysis 2'
        ]
        self.cases = [
            "case A",
            "case B",
            "case C",
            "case C qualitative",
            "case D",
            "case D not anhydride"
        ]


    def getPostData(self, calc, prop, method=None):
        return {'structure': "", 'melting_point': None}

    
    def makeDataRequest(self, url, structure, calc=None):
        _post = {'structure': structure}
        if self.melting_point != None:
            _post['melting_point'] = self.melting_point
        return self.request_logic(url, _post)

    
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
                    return self.results
                _retries += 1
            except Exception as e:
                logging.warning("Exception in calculator_epi.py: {}".format(e))
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
            logging.warning("epi server response status: {}".format(response.status_code))
            logging.warning("epi server response: {}".format(response.content))
            return False
        return True



    def get_mp_from_results(self, results):
        for data_obj in results['data']:
                if data_obj.get('prop') == 'melting_point':
                    logging.info("Found MP in EPI results..")
                    return float(data_obj['data'])
        return None


    def round_half_life(self, value):
        upper_bound = 1e3
        lower_bound = 1e-1

        if type(value) != float:
            value = float(value)

        if abs(value) > upper_bound or abs(value) < lower_bound:
            return "{:.2e}".format(value)
        else:
            return round(value, 2)


    def count_op_esters(self, child_nodes):
        """
        Returns number of products that have op ester 1 or 2 routes.
        """
        num_op_esters = 0
        for child in child_nodes:
            if child['route'] in self.op_esters:
                num_op_esters += 1
        return num_op_esters


    def handle_cases_a_and_b(self, response_obj, route, is_num_sites_1):
        """
        Handles OP Ester workflow cases A and B.
        """
        halfLifeValue = None

        if not is_num_sites_1:
            return None

        for data_obj in response_obj["data"]:
            # Accounting for cases A and B:
            if route == self.op_esters[0] and data_obj["prop"] == "Kb":
                halfLifeValue = self.round_half_life(data_obj["data"])
                break
            elif route == self.op_esters[1] and (data_obj["prop"] == "Ka" or data_obj["prop"] == "Kn"):
                halfLifeValue = self.round_half_life(data_obj["data"])
                break

        return halfLifeValue


    def handle_case_c(self, structure, response_obj, child_nodes, route, unique_schemes_count):
        """
        Handles OP Ester Workflow case C.
        """
        halfLifeValue = None
        case = None

        if not unique_schemes_count > 1:
            logging.warning("unique_schemes_count must be greater than 1 for case C. Returning None.")
            return None

        logging.info("unique_schemes_count > 1")
        logging.info("Case C: response_obj: {}".format(response_obj))

        if route in self.op_esters:

            logging.info("Check if OpProductCount > 4.")

            num_op_esters = self.count_op_esters(child_nodes)

            logging.info("OpProductCount: {}".format(num_op_esters))

            if num_op_esters > 4:
                # Return qualitative for all products, or just the OP Ester ones?
                logging.info("Return N/A (use qualitative values for HL).")
                case = "case C qualitative"
                halfLifeValue = None
                return halfLifeValue, case
            else:
                logging.info("Returning specific values if op ester 1 or 2.")
                if route == self.op_esters[0]:
                    logging.info("# op ester hydrolysis 1, returning Kb.")
                    kb_values = self.get_kb_values(response_obj)
                    logging.info("kb_values: {}".format(kb_values))
                    halfLifeValue = self.round_half_life(kb_values[0]['data'])
                    return halfLifeValue, None
                elif route == self.op_esters[1]:
                    logging.info("# op ester hydrolysis 2, returning Ka/n.")
                    ka_kn_values = self.get_ka_kn_values(response_obj)
                    halfLifeValue = self.round_half_life(ka_kn_values[0]['data'])
                    return halfLifeValue, None

        else:

            # if route == 'epoxide hydrolysis':
            if 'epoxide' in route.lower():
                logging.info("Sort Ka/n by atomNum.")

                sorted_response = self.sort_k_by_atom_number(response_obj, "Ka/n")

                # TODO: Validate sorted_response.

            elif route in self.cleaved_list:
                if route == 'anhydride hydrolysis' or route == 'carboxylic acid ester hydrolysis':
                    logging.info("Find functional group, assign to atomNum.")
                    func_groups = self.rdkit.get_functional_groups(route, structure)
                    logging.info("Functional groups for smiles: {}, route: {}: {}".format(structure, route, func_groups))
                    matched_data = self.match_atom_with_func_groups(response_obj, func_groups, ["Kb"])
                    if route == 'anhydride hydrolysis':
                        if len(matched_data) > 0 and "data" in matched_data[0]:
                            logging.info("Return Kb1, Kb2.")  # NOTE: think this happens on the epi wrapper side
                            halfLifeValue = self.round_half_life(matched_data[0]["data"])
                            logging.info("Setting half-life to {}.".format(halfLifeValue))
                    else:
                        # TODO: matched_data list gets sorted by atomNum, and (potentially, asking Lindsay) those values get
                        # split up amongst the child products.
                        logging.info("CAE - Sort Kb by atomNum.")
                        sorted_response = self.sort_k_by_atom_number(response_obj, "Kb")
                        logging.info("Sorted response: {}".format(sorted_response))

                else:
                    logging.info("Not CAE nor AH - Sort Kb by atomNum.")
                    sorted_response = self.sort_k_by_atom_number(response_obj, "Kb")
                    logging.info("Sorted response: {}".format(sorted_response))

            else:
                logging.info("No EH nor cleaved - Sort Kb by atomNum.")
                sorted_response = self.sort_k_by_atom_number(response_obj, "Kb")
                logging.info("Sorted response: {}".format(sorted_response))


        return halfLifeValue, case


    def handle_case_d(self, response_obj, route, unique_schemes_count):
        """
        Handles OP Ester Workflow case D.
        """
        halfLifeValue = None
        case = None

        if unique_schemes_count > 1:
            logging.warning("unique_schemes_count must be 1 or less for case D. Returning None.")
            return None

        for data_obj in response_obj["data"]:
            

            logging.info("\n\n")

            logging.info("unique_schemes_count <= 1")
            logging.info("Case D")

            logging.info("Half-life response data object: {}".format(data_obj))



            # Case D
            if "halogenated aliphatics" in route.lower():
                # return N/A
                logging.info("Case D halogenated aliphatics")
                halfLifeValue = "N/A"
                break
            elif "epoxide" in route.lower() and (data_obj["prop"] == "Ka" or data_obj["prop"] == "Kn"):
                # return Ka/n
                logging.info("Case D epoxide, returning Ka/n")
                halfLifeValue = self.round_half_life(data_obj["data"])
                break
            
            if route in self.cleaved_list:
                if "anhydride" in route.lower():
                    logging.info("Case D cleaved anhydride")
                    # return kb1 kb2 factored in (on epi wrapper side)
                    halfLifeValue = self.round_half_life(data_obj["data"])
                    break
                elif data_obj["prop"] == "Kb":
                    logging.info("case D not anhydride, returning Kb")

                    halfLifeValue1 = self.round_half_life(response_obj["data"][0]["data"])
                    halfLifeValue2 = self.round_half_life(response_obj["data"][1]["data"])
                    halfLifeValue = halfLifeValue1 + ", " + halfLifeValue2

                    case = "case D not anhydride"

                    break

            elif data_obj["prop"] == "Kb":
                logging.info("Case D not cleaved, returning Kb")
                # return Kb
                halfLifeValue = self.round_half_life(data_obj["data"])
                break

        return halfLifeValue, case


    def sort_k_by_atom_number(self, response_obj, prop):
        """
        Sorts Kb (or Ka/n) values by atomNum.
        For schemes: CAE, carbamates, epoxides, anhydrides
        """
        # Assuming sorting is ordering list of Kb values from half-life response to
        # be in ascending order by "atom_number" key.
        
        # So far, EPI seems to return them in order already?
        if not prop in ["Kb", "Ka/n"]:
            logging.error("sort_k_by_atom_number - Cannot filter prop '{}'. Must be Kb or Ka/n.")
            return None

        if not "data" in response_obj:
            logging.error("No 'data' key in half-life response object: {}".format(response_obj))
            return False

        # TODO: Create/filter list that's just Kb (or just Ka/n):
        filtered_hl_values = self.filter_by_property(response_obj, prop)

        sorted_hl_values = sorted(response_obj['data'], key=lambda x: int(x["atom_number"]))

        # What to do with the sorted values?

        return sorted_hl_values



    def filter_by_property(self, response_obj, prop):
        """
        Filters HL response by k value (Kb or Ka/n).
        """
        if not prop in ["Kb", "Ka/n"]:
            logging.error("filter_by_property - Cannot filter prop '{}'. Must be Kb or Ka/n.")
            return None

        if not "data" in response_obj:
            logging.error("Expected 'data' key in HL response.")
            return None

        logging.debug("Full half-life response: {}".format(response_obj))

        # filtered_response = [obj for obj in response_obj["data"] if obj.get("prop") == prop]
        filtered_response = []
        for hl_obj in response_obj["data"]:
            if prop == "Ka/n" and hl_obj.get("prop") in ["Ka", "Kn"]:
                filtered_response.append(hl_obj)
            elif prop == "Kb" and hl_obj.get("prop") == "Kb":
                filtered_response.append(hl_obj)

        logging.debug("Filtered half-life response: {}".format(filtered_response))

        return filtered_response


    def match_atom_with_func_groups(self, hl_response, func_group, props):
        """
        Matches atom number from HL response with functional group.
        """
        if not "data" in hl_response:
            logging.warning("Data key not in hl_response.")
            return False

        matched_data = []

        logging.info("HL Response: {}.".format(hl_response))

        for data_obj in hl_response["data"]:
            
            atom_number = data_obj["atom_number"]
            logging.info("Atom number: {}.".format(atom_number))

            if not atom_number:
                logging.info("Atom number is null, skipping to next one.")
                continue

            for group_number in func_group:

                logging.info("Group number: {}.".format(group_number))

                # if int(group_number) == int(atom_number) and data_obj["prop"] == prop:
                if int(group_number) == int(atom_number) and data_obj["prop"] in props:
                    logging.info("Matched atom with group number.")
                    matched_data.append(data_obj)

        logging.info("Matched data: {}.".format(matched_data))

        if len(matched_data) > 0:
            logging.warning("No matches for functional group and atom numbers. Returning blank HL object.")
            return [{"data": None}]  # returns empty HL response

        return matched_data


    def get_kb_values(self, response_obj):
        """
        Returns data for Kb values.
        """
        kb_data = []
        for data_obj in response_obj.get("data"):
            if data_obj["prop"] == "Kb":
                kb_data.append(data_obj)
        return kb_data


    def get_ka_kn_values(self, response_obj):
        """
        Returns data for Kb values.
        """
        ka_kn_data = []
        for data_obj in response_obj.get("data"):
            if data_obj["prop"] == "Ka" or data_obj["prop"] == "Kn":
                ka_kn_data.append(data_obj)
        return ka_kn_data


    def assign_qualitative_values(self, child_nodes):
        """
        Using qualitative descriptor for half life values.
        """
        qsar_responses = []
        for child_obj in child_nodes:
            qsar_response = {}
            qsar_response.update(child_obj)
            qsar_response["data"] = None
            qsar_response["prop"] = "qsar"
            qsar_response["valid"] = False
            qsar_responses.append(qsar_response)
        return qsar_responses


    def assign_both_kb(self, child_nodes, halfLifeValue):
        """
        Two Kb values are returned and the first half of children get the first,
        and the second half get the second half-life.
        """

        
        # NOTE: This is assuming that two Kb values are returned (but what about > 2 Kb values?)


        halfLifeValue1 = halfLifeValue.split(",")[0]
        halfLifeValue2 = halfLifeValue.split(",")[1]
        i = 0
        qsar_responses = []
        for child_obj in child_nodes:
            qsar_response = {}
            qsar_response.update(child_obj)
            qsar_response["prop"] = "qsar"
            qsar_response["valid"] = True

            # qsar_responses.append(qsar_response)
            if i < len(child_nodes) / 2:
                # assign first value
                qsar_response["data"] = halfLifeValue1
            else:
                # assign second value
                qsar_response["data"] = halfLifeValue2

            qsar_responses.append(qsar_response)

            i += 1

        return qsar_responses



    def determine_halflife(self, structure, response_obj, child_nodes, route, is_num_sites_1, unique_schemes_count):
        """
        Determines halflife based on response.
        """
        halfLifeValue = None
        case = None

        if is_num_sites_1:
            # Cases A and B
            halfLifeValue = self.handle_cases_a_and_b(response_obj, route, is_num_sites_1)
        else:
            if unique_schemes_count > 1:
                # Case C
                halfLifeValue, case = self.handle_case_c(structure, response_obj, child_nodes, route, unique_schemes_count)
            else:
                # Case D
                halfLifeValue, case = self.handle_case_d(response_obj, route, unique_schemes_count)

        if not halfLifeValue:
            logging.warning("halfLifeValue not set, getting data from response_obj: {}".format(response_obj))
            halfLifeValue = self.round_half_life(response_obj["data"][0]["data"])

        return halfLifeValue, case


    def make_qsar_request(self, request_dict):
        """
        Makes requests to epi suite for half-lives.
        """

        structure = request_dict.get("filtered_smiles")
        unique_schemes_count = int(request_dict.get("uniqueSchemesCount"))  
        product_count = int(request_dict.get("productCount"))
        child_nodes = request_dict.get("childNodes")  # children of a single/given parent
        qsar_responses = []

        logging.info("\n\n\n Starting main QSAR request loop.")

        for child_obj in child_nodes:

            logging.info("~Child object: {}".format(child_obj))

            route = child_obj.get("routes").lower()
            route_url = None
            num_sites = None
            is_num_sites_1 = False  # whether num_sites is == 1 or > 1

            if not route in list(self.qsar_request_map.keys()):
                # raise Exception("Route not found: {}.".format(route))
                logging.warning("Route not found: {}".format(route))
                qsar_response = {}
                qsar_response.update(child_obj)
                qsar_response["data"] = None
                qsar_response["prop"] = "qsar"
                qsar_response["valid"] = False
                qsar_responses.append(qsar_response)
                continue

            route_url = self.qsar_request_map[route]
            url = self.baseUrl.replace("estimated", "") + route_url

            if not route in self.cleaved_list:
                logging.info("Route not in cleaved list. Product count: {}".format(product_count))
                if product_count > 1:
                    is_num_sites_1 = False
                else:
                    is_num_sites_1 = True

            elif route in self.cleaved_list and not route in self.op_esters:
                logging.info("Route in cleaved list but not OP Ester. Product count: {}".format(product_count))
                if product_count > 2:
                    is_num_sites_1 = False
                else:
                    is_num_sites_1 = True

            elif route in self.cleaved_list and route in self.op_esters:
                logging.info("Route in cleaved list and an OP Ester. Product count: {}".format(product_count))
                if product_count > 4:
                    logging.info("Product count > 4.")
                    is_num_sites_1 = False
                    # qsar_response = {}
                    # qsar_response.update(child_obj)
                    # qsar_response["data"] = None
                    # qsar_response["prop"] = "qsar"
                    # qsar_response["valid"] = True
                    # qsar_responses.append(qsar_response)
                    # continue
                else:
                    logging.info("Product count <= 4.")
                    is_num_sites_1 = True

            else:
                is_num_sites_1 = True

            logging.info(">>> Is number of sites one?: {}".format(is_num_sites_1))
            logging.info("Request to EPI for half life:\nURL:{}\nStructure:{}".format(url, structure))

            response = requests.post(url, data=json.dumps({'structure': structure}), headers=self.headers)

            if response.status_code != 200:
                logging.warning("Error requesting half-life data from EPI Suite.\nStatus code: {}\nContent: {}".format(response.status_code, response.content))
                qsar_response = {}
                qsar_response.update(child_obj)
                qsar_response["data"] = None
                qsar_response["prop"] = "qsar"
                qsar_response["valid"] = False
                qsar_responses.append(qsar_response)
                continue

            # try:
            response_obj = json.loads(response.content)

            if not response_obj.get("data") or len(response_obj.get("data")) < 1:
                raise Exception("Half-life response does not have excepted 'data' key or size is unexpected.\nResponse: {}".format(response_obj))
            
            # NOTE: If OP Ester 1/2 and num_sites <= 1, pick specific half-life value from set.
            # Get Kb if OP Ester 1, and Ka/n for OP Ester 2

            halfLifeValue, case = self.determine_halflife(structure, response_obj, child_nodes, route, is_num_sites_1, unique_schemes_count)

            if case == "case D not anhydride":
                # Assign HL1 to first half of children, HL2 to second half, no need to 
                # continue looping and requesting.
                logging.info("case D not anhydride.")
                qsar_responses = self.assign_both_kb(child_nodes, halfLifeValue)
                return qsar_responses
            # elif case == "case C qualitative":
            #     # Uses qualitative descriptor for half life values
            #     logging.info("case C qualitative.")
            #     qsar_responses = self.assign_qualitative_values(child_nodes)
            #     return qsar_responses

            qsar_response = {}
            qsar_response.update(child_obj)
            qsar_response["data"] = halfLifeValue
            qsar_response["prop"] = "qsar"
            qsar_response["valid"] = True

            logging.info("QSAR response: {}\n\n".format(qsar_response))

            qsar_responses.append(qsar_response)

            # except Exception as e:
            #     logging.warning("Error parsing QSAR response: {}".format(e))
            #     qsar_response = {}
            #     qsar_response["error"] = "Error getting QSAR data from EPI Suite"
            #     qsar_response["valid"] = False
            #     qsar_responses.append(qsar_response)


        return qsar_responses


    def data_request_handler(self, request_dict):
        """
        Makes requests to the EPI Suite server
        """
        
        _filtered_smiles = ''
        _response_dict = {}

        # fill any overlapping keys from request:
        for key in request_dict.keys():
            if not key == 'nodes':
                _response_dict[key] = request_dict.get(key)
        _response_dict.update({'request_post': request_dict, 'method': None})

        try:
            _filtered_smiles = SMILESFilter().parseSmilesByCalculator(request_dict['chemical'], request_dict['calc']) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            _response_dict.update({
                'data': "Cannot filter SMILES",
                'valid': False
            })
            return _response_dict

        # Handle QSAR request or continue to the usual p-chem stuff
        if request_dict.get('prop') == 'qsar':
            
            request_dict['filtered_smiles'] = _filtered_smiles

            _result_obj = self.make_qsar_request(request_dict)

            # TODO: Account for not valid result_obj

            _response_dict['data'] = _result_obj

            # _response_dict.update(_result_obj['data'])
            _response_dict['valid'] = True

            return _response_dict

        try:

            _get_mp = request_dict.get('prop') == 'water_sol' or request_dict.get('prop') == 'vapor_press'
            
            if _get_mp:
                self.melting_point = self.get_melting_point(_filtered_smiles, 
                                        request_dict.get('sessionid'), self)
            else:
                self.melting_point = None

            _result_obj = self.makeDataRequest(self.baseUrl, _filtered_smiles, request_dict['calc']) # make call for data!

            if _get_mp and not self.melting_point:
                # MP not found from measured or test, getting from results,
                # and requesting data again with set MP..
                self.melting_point = self.get_mp_from_results(_result_obj)
                _result_obj = self.makeDataRequest(self.baseUrl, _filtered_smiles, request_dict['calc'])  # Make request using MP

            _response_dict.update(_result_obj)
            _response_dict['valid'] = True
        
            return _response_dict

        except Exception as err:
            logging.warning("Exception occurred getting {} data: {}".format(err, request_dict['calc']))
            _response_dict.update({
                'data': "Cannot reach EPI calculator",
                'valid': False
            })
            return _response_dict
