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

    def flatten(self, l):
        """
        Flattens a list of lists.
        """
        return [item for sublist in l for item in sublist]

    def increment_atom_number(self, atom_list):
        """
        GetSubstructMatches uses 0-start index hydrowin use 1-start index.
        this functions adds +1 to all list items in the return atom number list in the FindFG function to match better
        """
        updated_atoms = []
        for atom in atom_list:
            new_atom = atom + 1
            updated_atoms.append(new_atom)
        return updated_atoms

    def get_functional_groups_anhydride(self, smiles):
        """
        List of atoms in anhydride functional group.
        """
        mol = Chem.MolFromSmiles(smiles)
        anhydride_atom = self.flatten(list(Chem.Mol.GetSubstructMatches(mol, self.CarbAnhydride, uniquify=True)))
        logging.info("Anhydride group: {}".format(anhydride_atom))
        anhydride_atom = self.increment_atom_number(anhydride_atom)
        logging.info("Updated Anhydride group: {}".format(anhydride_atom))
        return anhydride_atom

    def get_functional_groups_cae(self, smiles):
        """
        List of atoms in carboxylic acid ester functional group.
        """
        mol = Chem.MolFromSmiles(smiles)
        CAE_atom = self.flatten(list(Chem.Mol.GetSubstructMatches(mol, self.CAE, uniquify=True)))
        logging.info("CAE group: {}".format(CAE_atom))
        CAE_atom = self.increment_atom_number(CAE_atom)
        logging.info("Updated CAE group: {}".format(CAE_atom))
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
        logging.info("Counting op esters, child_nodes: {}".format(child_nodes))
        num_op_esters = 0
        for child in child_nodes:
            if child.get("routes").lower() in self.op_esters:
                num_op_esters += 1
        return num_op_esters


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


    def match_atom_with_func_groups(self, route, structure, hl_response, props):
        """
        Matches atom number from HL response with functional group.
        """
        if not "data" in hl_response:
            logging.warning("Data key not in hl_response.")
            return False

        func_group = self.rdkit.get_functional_groups(route, structure)

        matched_data = []

        logging.warning("HL Response: {}.".format(hl_response))
        logging.warning("Functional groups: {}".format(func_group))

        for data_obj in hl_response["data"]:
            
            atom_number = data_obj["atom_number"]
            logging.warning("Atom number: {}.".format(atom_number))

            if not atom_number:
                logging.warning("Atom number is null, skipping to next one.")
                continue

            for group_number in func_group:

                logging.warning("Group number: {}.".format(group_number))

                # if int(group_number) == int(atom_number) and data_obj["prop"] == prop:
                if int(group_number) == int(atom_number) and data_obj["prop"] in props:
                    logging.warning("Matched atom with group number.")
                    matched_data.append(data_obj)

        logging.info("Matched data: {}.".format(matched_data))

        if len(matched_data) < 1:
            logging.warning("No matches for functional group and atom numbers. Returning blank HL object.")
            return []  # returns empty HL response

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
        for child_obj in child_nodes:
            child_obj["data"] = None
            child_obj["prop"] = "qsar"
            child_obj["valid"] = False
        return child_nodes


    def is_num_sites_1(self, child_obj, product_count, route):
        """
        Returns boolean if number of sites is 1 or not.
        """
        is_one = False  # whether num_sites is == 1 or > 1

        if not route in self.cleaved_list:
            logging.info("Route not in cleaved list. Product count: {}".format(product_count))
            if product_count > 1:
                is_one = False
            else:
                is_one = True
        elif route in self.cleaved_list and not route in self.op_esters:
            logging.info("Route in cleaved list but not OP Ester. Product count: {}".format(product_count))
            if product_count > 2:
                is_one = False
            else:
                is_one = True
        elif route in self.cleaved_list and route in self.op_esters:
            logging.info("Route in cleaved list and an OP Ester. Product count: {}".format(product_count))
            if product_count > 4:
                logging.info("Product count > 4.")
                is_one = False
            else:
                logging.info("Product count <= 4.")
                is_one = True
        else:
            is_one = True

        return is_one


    def is_op_ester(self, route):
        """
        Returns boolean indicating if route is an op ester or not.
        """
        if route.lower() in self.op_esters:
            return True
        else:
            return False


    def sort_products_by_case(self, parent, unique_schemes_count, product_count, child_nodes):
        """
        Loops child nodes and organizes them by case before determining HL.
        """
        qsar_map = {
            ""
        }
        for child_obj in child_nodes:

            logging.info("ORIGINAL CHILD OBJ: {}".format(child_obj))

            route = child_obj.get("routes").lower()
            is_one = self.is_num_sites_1(child_obj, product_count, route)
            op_ester = self.is_op_ester(route)

            logging.info("Route: {}".format(route))
            logging.info("Is one: {}".format(is_one))
            logging.info("Is OP ester: {}".format(op_ester))
            logging.info("Unique schemes: {}".format(unique_schemes_count))

            if is_one:
                if op_ester:
                    # case B
                    child_obj["case"] = "B"
                    child_obj["path"] = "1"  # NOTE: case B only has one path: single-site OP Esters
                else:
                    # NOTE: case A (no op esters, right? that's case B)
                    child_obj["case"] = "A"
                    child_obj["path"] = self.determine_path(child_obj, route, child_nodes)
            else:
                if unique_schemes_count > 1:
                    # case C
                    child_obj["case"] = "C"
                    child_obj["path"] = self.determine_path(child_obj, route, child_nodes)
                else:
                    # case D
                    child_obj["case"] = "D"
                    child_obj["path"] = self.determine_path(child_obj, route, child_nodes)

            logging.info("UPDATED CHILD OBJ: {}".format(child_obj))

        return child_nodes


    def determine_path(self, child_obj, route, child_nodes):
        """
        Determine path for a product's case.
        """
        if not child_obj.get("case"):
            logging.warning("determine_path() - 'case' not in child_obj.")
            return False

        path = None

        if child_obj["case"] == "A":
            path = self.handle_case_a_path(route)
        elif child_obj["case"] == "B":
            path = "1"  # NOTE: only one path for case B
        elif child_obj["case"] == "C":
            path = self.handle_case_c_path(route, child_nodes)
        elif child_obj["case"] == "D":
            path = self.handle_case_d_path(route)

        return path


    def handle_case_a_path(self, route):
        """
        Case A - for single-site non-op esters.
        """
        path = None
        if route == "epoxide":
            path = "1"
        elif route in self.cleaved_list:
            if "anhydride" in route:
                path = "4"
            else:
                path = "5"
        elif "halogenated aliphatics" in route:
            path = "2"
        else:
            path = "3"
        return path


    def handle_case_c_path(self, route, child_nodes):
        """
        Case C - for multi-site with unique schemes.
        """
        path = None
        if route in self.op_esters:
            num_op_esters = self.count_op_esters(child_nodes)
            if num_op_esters > 4:
                path = "1"
            else:
                path = "2"
        elif "epoxide" in route:
            path = "3"
        elif route in self.cleaved_list:
            path = "5"
        else:
            path = "4"
        return path

    def handle_case_d_path(self, route):
        """
        Case D - for multi-site with single scheme.
        """
        path = None
        if "halogenated aliphatics" in route:
            path = "1"
        elif "epoxide" in route:
            path = "2"
        elif route in self.cleaved_list:
            path = "3"
        else:
            path = "4"
        return path


    def group_products(self, child_nodes):
        """
        Loops child products and groups them by their case and path/scheme.
        """
        # Creating key for group that's case + path, e.g., "A2" = case A path 2
        grouped_products = {}
        for child_obj in child_nodes:
            key = child_obj["case"] + child_obj["path"]
            if key in list(grouped_products.keys()):
                grouped_products[key].append(child_obj)
            else:
                grouped_products[key] = [child_obj]
        return grouped_products


    def validate_path_routes(self, child_obj_list):
        """
        Validates that all the child_obj's for a given path (e.g., "A5") 
        are the same route.
        """
        if len(child_obj_list) < 1:
            logging.warning("Cannot check routes, child_obj_list less than one.")
            return False

        init_route = child_obj_list[0].get("routes")
        
        for child_obj in child_obj_list[1:]:
            route = child_obj.get("routes")

            if self.is_op_ester(init_route) and not self.is_op_ester(route):
                logging.warning("OP ester routes of the same path do not match.")
                return False

            elif not self.is_op_ester(init_route) and route != init_route:
                logging.warning("Routes of the same path do not match.")
                return False

        return True


    def get_qsar_for_products(self, parent, grouped_products):
        """
        Makes request to EPI for QSAR data.

        grouped_products example:
            {
                "A5": [{child_obj}, {child_obj}, ..],
                "B1": [{child_obj}, {child_obj}, ..]
            }

        NOTE: Each key indicating a path (e.g., "A5") should have the same route/scheme.
        """
        logging.info("get_qsar_for_products() grouped_products: {}".format(grouped_products))

        all_products_list = []

        for path_key, child_obj_list in grouped_products.items():

            logging.info("Path key: {}\nchild_obj_list: {}".format(path_key, child_obj_list))

            case = path_key[0]
            path = path_key[1]
            route = child_obj_list[0]["routes"].lower()  # routes for child_obj_list should all be the same
            route_endpoint = self.qsar_request_map[route]
            url = self.baseUrl.replace("estimated", "") + route_endpoint

            logging.info("Path key: {}\nRoute: {}\nUrl: {}".format(path_key, route, url))

            response_obj = {
                "status": False,
                "qsar_response": None
            }

            if path_key in ["A2", "C1", "D1"]:
                logging.info("Skipping request for case: {}, path: {}, assigning qualitative values.".format(case, path))
                child_obj_list = self.assign_qualitative_values(child_obj_list)
                all_products_list += child_obj_list
                continue

            try:

                logging.info("Making QSAR request to EPI.")

                response = requests.post(url, data=json.dumps({'structure': parent}), headers=self.headers)

                if response.status_code != 200:
                    logging.warning("Error requesting half-life data from EPI Suite.\nStatus code: {}\nContent: {}".format(response.status_code, response.content))
                    for child_obj in child_obj_list:
                        child_obj["error"] = "Error requesting half-life data from EPI."
                        child_obj["prop"] = "qsar"
                        child_obj["valid"] = False
                    all_products_list += child_obj_list
                    continue

                response_obj = json.loads(response.content)

                if not response_obj.get("data") or len(response_obj.get("data")) < 1:
                    logging.warning("Error parsing half-life data from EPI response.")
                    for child_obj in child_obj_list:
                        child_obj["error"] = "Error parsing half-life data from EPI response."
                        child_obj["prop"] = "qsar"
                        child_obj["valid"] = False
                    all_products_list += child_obj_list
                    continue
                
                # Assigns data to products in list based on case and path.
                child_obj_list = self.handle_hl_response(response_obj, parent, route, case, path, child_obj_list)

                all_products_list += child_obj_list

            except Exception as e:
                logging.warning("Error making QSAR request: {}".format(e))
                for child_obj in child_obj_list:
                    # child_obj["data"] = None
                    child_obj["error"] = "Error making request to EPI for half-life."
                    child_obj["prop"] = "qsar"
                    child_obj["valid"] = False
                continue

        return all_products_list


    def handle_hl_response(self, response_obj, parent, route, case, path, child_obj_list):
        """
        Assigns HL to product based on case and path.
        """
        logging.info("HL Response: {}".format(response_obj))

        site_type = ""

        if case in ["A", "B"]:
            site_type = "single"
        elif case in ["C", "D"]:
            site_type = "multi"

        logging.info("Site type: {}".format(site_type))

        num_hls = len(response_obj["data"])

        logging.info("Number of HLs: {}".format(num_hls))
            
        if case == "A":
            if path == "1":
                return self.hl_result_pattern("Ka/n", child_obj_list, response_obj)
            elif path == "2":
                return self.assign_qualitative_values(child_obj_list)
            elif path == "4":
                return self.handle_functional_group_case(route, parent, response_obj, child_obj_list)
            else:
                return self.hl_result_pattern("Kb", child_obj_list, response_obj)
        elif case == "B":
            return self.handle_op_ester_values(response_obj, child_obj_list)
        elif case == "D":
            child_obj_list = self.handle_case_d_results(path, response_obj, child_obj_list)
        elif case == "C":
            child_obj_list = self.handle_case_c_results(parent, path, route, response_obj, child_obj_list)

        return child_obj_list


    def handle_op_ester_values(self, response_obj, child_obj_list):
        """
        Assigns HL values based on OP ester 1 or 2.
        """
        op1_hl = None
        op2_hl = None

        # Gets HL values for OP ester 1 and 2:
        for data_obj in response_obj.get("data"):
            if data_obj["prop"] == "Kb":
                op1_hl = data_obj.get("data")
            elif data_obj["prop"] == "Ka" or data_obj["prop"] == "Kn":
                op2_hl = data_obj.get("data")

        # Splits up op ester 1 and 2 values:
        for child_obj in child_obj_list:
            route = child_obj.get("routes").lower()
            if route == self.op_esters[0]:
                child_obj["data"] = self.round_half_life(op1_hl)
            elif route == self.op_esters[1]:
                child_obj["data"] = self.round_half_life(op2_hl)

        return child_obj_list


    def handle_functional_group_case(self, route, parent, response_obj, child_obj_list):
        """
        Finds HLs for any atom numbers that match ones returned 
        by the functional groups. If number of sites is one, returns
        an HL value.
        """
        matched_data = self.match_atom_with_func_groups(route, parent, response_obj, ["Kb"])
        if len(matched_data) == 1 and "data" in matched_data[0]:
            logging.warning("Matched Data: {}".format(matched_data))
            return self.hl_result_pattern("Kb", child_obj_list, {"data": matched_data})
        logging.warning("Matched data not equal to one: {}. Using qualitative values.".format(matched_data))
        return self.assign_qualitative_values(child_obj_list)


    def handle_case_d_results(self, path, response_obj, child_obj_list):
        """
        Assigns HL values to products based on case D paths.
        """
        logging.info("Case: D{}".format(path))
        if path == "1":
            return self.assign_qualitative_values(child_obj_list)
        elif path == "2":
            return self.hl_result_pattern("Ka/n", child_obj_list, response_obj)
        elif path in ["3", "4"]:
            return self.hl_result_pattern("Kb", child_obj_list, response_obj)


    def handle_case_c_results(self, parent, path, route, response_obj, child_obj_list):
        """
        Assigns HL values to products based on case C paths.
        """
        logging.info("Case: C{}".format(path))
        if path == "1":
            return self.assign_qualitative_values(child_obj_list)
        elif path == "2":
            return self.handle_op_ester_values(response_obj, child_obj_list)
        elif path == "3":
            return self.hl_result_pattern("Ka/n", child_obj_list, response_obj)
        elif path == "4" or path == "5":
            return self.hl_result_pattern("Kb", child_obj_list, response_obj)


    def hl_result_pattern(self, sort_prop, child_obj_list, response_obj):
        """
        Consolidating similar HL value assignment logic that accounts
        for number of HL values. Shared across various cases.
        """
        sorted_response = self.sort_k_by_atom_number(response_obj, sort_prop)
        if len(sorted_response) == 1:
            for child_obj in child_obj_list:
                child_obj["data"] = self.round_half_life(sorted_response[0]["data"])
            return child_obj_list
        elif len(sorted_response) > 1:
            return self.split_hl_values(child_obj_list, sorted_response, "Kb")
        else:
            raise Exception("I'm not even supposed to be here.")


    def split_hl_values(self, child_obj_list, sorted_response, sort_prop="Kb"):
        """
        Splits up HL values across products.

        NOTE: This assumes two HLs and an even number of products!
        """
        # sorted_response = self.sort_k_by_atom_number(response_obj, sort_prop)  # NOTE: Now called in hl_result_pattern()
        logging.info("Sorted response: {}".format(sorted_response))

        num_products = len(child_obj_list)
        num_hls = len(sorted_response)
        mid_point = int(num_products / 2)

        logging.info("Mid point: {}\nMid point type: {}".format(mid_point, type(mid_point)))

        logging.info("Number of products: {}\nNumber of HLs: {}".format(num_products, num_hls))

        # NOTE: This setup assumes that there are only 2 HLs

        for child_obj in child_obj_list[0:mid_point]:
            child_obj["data"] = self.round_half_life(sorted_response[0]["data"])

        for child_obj in child_obj_list[mid_point:]:
            child_obj["data"] = self.round_half_life(sorted_response[1]["data"])

        return child_obj_list


    def make_qsar_request(self, request_dict):
        """
        Makes requests to epi suite for half-lives.
        """

        # structure = request_dict.get("filtered_smiles")
        parent = request_dict.get("filtered_smiles")
        unique_schemes_count = int(request_dict.get("uniqueSchemesCount"))  
        product_count = int(request_dict.get("productCount"))
        child_nodes = request_dict.get("childNodes")  # children of a single/given parent
        qsar_responses = []

        logging.info("\n\nStarting main QSAR request loop.\n\n")

        child_nodes = self.sort_products_by_case(parent, unique_schemes_count, product_count, child_nodes)

        grouped_products = self.group_products(child_nodes)

        logging.info("GROUPED PRODUCTS: {}".format(grouped_products))

        qsar_responses = self.get_qsar_for_products(parent, grouped_products)

        logging.info("UPDATED CHILD NODES TO USE TO DETERMINE HL REQUESTS AND VALUE ASSIGNMENT: {}".format(child_nodes))

        logging.info("QSAR RESPONSES: {}".format(qsar_responses))

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
