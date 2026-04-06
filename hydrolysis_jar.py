import logging
import json
import requests
import os
import re
import math

from .calculator_rdkit import RdkitCalc


# ---------------------------------------------------------------------------
# Half-life conversion helpers (ported from EPIReader.cs)
# ---------------------------------------------------------------------------

def _rate_to_half_life_days(rate: float, is_neutral: bool = False) -> float:
    """
    Convert a hydrolysis rate constant to a half-life in days.

    For Kb (base-catalyzed):  t½ = 0.6931 / (k * 1e-7) / 86400
    For Kn (neutral):         t½ = 0.6931 / k / 86400

    NOTE: The C# original applies the 1e-7 factor for Kb but not Kn — mirrored here.
    """
    try:
        if is_neutral:
            return 0.6931 / rate / 86400.0
        else:
            return 0.6931 / (rate * 1.0e-7) / 86400.0
    except (ZeroDivisionError, ValueError):
        return math.nan


def _parse_rate_value(token: str) -> float:
    """
    Safely parse a scientific notation string to float.
    """
    try:
        return float(token)
    except (ValueError, TypeError):
        return math.nan



def parse_hydrolysis_output(output_text: str, prop_name: str, acid_base: str = "Kb") -> list[dict]:
    """
    Searches the EPI output text for a functional group label (prop_name, e.g.
    'CARBAMATE', 'ESTER') then finds every 'Kb hydrolysis at atom # N' line
    beneath it and converts the rate constant to a half-life in days.

    Returns a list of dicts matching the existing response_obj["data"] shape:
        [{"prop": "Kb", "data": "<days>", "units": "days", "atom_number": "3"}, ...]
    """
    results = []

    # Find the starting position of the functional group label (case-insensitive)
    search_start = output_text.lower().find(prop_name.lower())
    if search_start < 0:
        return results

    search_str = f"{acid_base} hydrolysis at "
    idx = search_start

    while True:
        idx = output_text.lower().find(search_str.lower(), idx)
        if idx < 0:
            break

        line_end = output_text.find("\n", idx)
        line = output_text[idx:line_end].strip()

        # Line format: "Kb hydrolysis at atom #  3:  4.680E+000  L/mol-sec"
        parts = line.split(":")
        if len(parts) < 2:
            idx += len(search_str)
            continue

        # Extract atom number from "Kb hydrolysis at atom #  3"
        atom_part = parts[0]
        atom_num = ""
        if "#" in atom_part:
            atom_num = atom_part.split("#")[-1].strip()

        # Extract rate value — first whitespace-delimited token after the colon
        value_tokens = parts[1].strip().split()
        if not value_tokens:
            idx += len(search_str)
            continue

        rate = _parse_rate_value(value_tokens[0])
        half_life = _rate_to_half_life_days(rate, is_neutral=False)

        results.append({
            "prop": acid_base,
            "data": str(half_life) if not math.isnan(half_life) else None,
            "units": "days",
            "atom_number": atom_num,
        })

        idx += len(search_str)

    return results


def parse_phosphate_thiophosphate_output(output_text: str) -> list[dict]:
    """
    Handles both Kn (neutral) and Kb (base-catalyzed) rate constants from
    phosphate/thiophosphate EPI output. Applies different math for each.

    Line format:
        Kn = 1.516e-007/sec = 9.097e-006/min
        Kb = 0.008675/M-sec = 0.5205/M-min
    """
    results = []

    for line in output_text.splitlines():
        stripped = line.strip()
        for prop_type in ("Kn", "Kb"):
            if not re.match(r"^" + prop_type, stripped):
                continue

            # Split on '=' and '/' to isolate the /sec value
            # "Kn = 1.516e-007/sec = ..." → tokens[1] is "1.516e-007"
            tokens = re.split(r"[=/]", stripped)
            if len(tokens) < 2:
                continue

            rate = _parse_rate_value(tokens[1].strip())
            is_neutral = (prop_type == "Kn")
            half_life = _rate_to_half_life_days(rate, is_neutral=is_neutral)

            results.append({
                "prop": prop_type,
                "data": str(half_life) if not math.isnan(half_life) else None,
                "units": "days",
                "atom_number": None,
            })

    return results


def parse_anhydride_output(output_text: str) -> list[dict]:
    """
    Looks for the 'Total Kb for pH' line which gives the combined rate constant
    for both ester sites in an anhydride, then converts to half-life in days.

    Line format:
        Total Kb for pH > 8 at 25 deg C :  2.570E+003  L/mol-sec
    """
    results = []

    for line in output_text.splitlines():
        stripped = line.strip()
        if not re.match(r"^Total Kb for pH", stripped, re.IGNORECASE):
            continue

        tokens = stripped.split(":")
        if len(tokens) < 2:
            break

        value_tokens = tokens[1].strip().split()
        if not value_tokens:
            break

        rate = _parse_rate_value(value_tokens[0])
        half_life = _rate_to_half_life_days(rate, is_neutral=False)

        results.append({
            "prop": "Kb",
            "data": str(half_life) if not math.isnan(half_life) else None,
            "units": "days",
            "atom_number": None,
        })
        break  # Only one Total Kb line expected

    return results


def parse_output_by_route(output_text, route):
    """
    Chooses the correct parser based on the hydrolysis route.
    """
    route = route.lower()

    if "phosphate" in route or "organophosphorus" in route:
        return parse_phosphate_thiophosphate_output(output_text)

    elif "anhydride" in route:
        return parse_anhydride_output(output_text)

    elif "carbamate" in route:
        return parse_hydrolysis_output(output_text, prop_name="CARBAMATE", acid_base="Kb")

    elif "ester" in route:
        return parse_hydrolysis_output(output_text, prop_name="ESTER", acid_base="Kb")

    elif "epoxide" in route:
        return parse_hydrolysis_output(output_text, prop_name="EPOXIDE", acid_base="Ka")

    elif "alkylhalide" in route or "halogenated" in route:
        # Alkyl halides use Kn (neutral) hydrolysis
        return parse_hydrolysis_output(output_text, prop_name="ALKYL HALIDE", acid_base="Kn")

    else:
        logging.warning(f"parse_output_by_route: unrecognized route '{route}', returning empty.")
        return []



class Hydrolysis:

    def __init__(self):

        self.rdkit = RdkitCalc()

        self.baseUrl = os.environ['CTS_EPI_SERVER']

        self.headers = {'Content-Type': 'application/json'}

        self.qsar_request_map = {
            'halogenated aliphatics: elimination': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (no adjacent x)': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (vicinal x)': 'hydrolysis/alkylhalide',
            'halogenated aliphatics: nucleophilic substitution (geminal x)': 'hydrolysis/alkylhalide',
            'epoxide hydrolysis': 'hydrolysis/epoxide',
            'organophosphorus ester hydrolysis 1': 'hydrolysis/phosphate',
            'organophosphorus ester hydrolysis 2': 'hydrolysis/phosphate',
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

    def round_half_life(self, value):
        upper_bound = 1e3
        lower_bound = 1e-1

        if not value:
            return 0

        if type(value) != float:
            value = float(value)

        if abs(value) > upper_bound or abs(value) < lower_bound:
            return "{:.2e}".format(value)
        else:
            return round(value, 2)

    def count_op_esters(self, child_nodes):
        logging.info("Counting op esters, child_nodes: {}".format(child_nodes))
        num_op_esters = 0
        for child in child_nodes:
            if child.get("routes").lower() in self.op_esters:
                num_op_esters += 1
        return num_op_esters

    def sort_k_by_atom_number(self, response_obj, prop):
        if not prop in ["Kb", "Ka/n"]:
            logging.error("sort_k_by_atom_number - Cannot filter prop '{}'.".format(prop))
            return None

        if not "data" in response_obj:
            logging.error("No 'data' key in half-life response object: {}".format(response_obj))
            return False

        filtered_hl_values = self.filter_by_property(response_obj, prop)

        def sorting_func(item):
            if "atom_number" in item and item["atom_number"] is not None:
                return int(item["atom_number"])
            return 0

        sorted_hl_values = sorted(response_obj['data'], key=sorting_func)
        return sorted_hl_values

    def filter_by_property(self, response_obj, prop):
        if not prop in ["Kb", "Ka/n"]:
            logging.error("filter_by_property - Cannot filter prop '{}'.".format(prop))
            return None

        if not "data" in response_obj:
            logging.error("Expected 'data' key in HL response.")
            return None

        filtered_response = []
        for hl_obj in response_obj["data"]:
            if prop == "Ka/n" and hl_obj.get("prop") in ["Ka", "Kn"]:
                filtered_response.append(hl_obj)
            elif prop == "Kb" and hl_obj.get("prop") == "Kb":
                filtered_response.append(hl_obj)

        return filtered_response

    def assign_qualitative_values(self, child_nodes):
        for child_obj in child_nodes:
            child_obj["data"] = None
            child_obj["prop"] = "qsar"
            child_obj["valid"] = False
        return child_nodes

    def is_num_sites_1(self, child_obj, product_count, route):
        is_one = False

        if not route in self.cleaved_list:
            is_one = product_count <= 1
        elif route in self.cleaved_list and route not in self.op_esters:
            is_one = product_count <= 2
        elif route in self.cleaved_list and route in self.op_esters:
            is_one = product_count <= 4
        else:
            is_one = True

        return is_one

    def is_op_ester(self, route):
        return route.lower() in self.op_esters

    def sort_products_by_case(self, parent, unique_schemes_count, product_count, child_nodes):
        for child_obj in child_nodes:
            route = child_obj.get("routes").lower()
            is_one = self.is_num_sites_1(child_obj, product_count, route)
            op_ester = self.is_op_ester(route)

            if is_one:
                if op_ester:
                    child_obj["case"] = "B"
                    child_obj["path"] = "1"
                else:
                    child_obj["case"] = "A"
                    child_obj["path"] = self.determine_path(child_obj, route, child_nodes)
            else:
                if unique_schemes_count > 1:
                    child_obj["case"] = "C"
                    child_obj["path"] = self.determine_path(child_obj, route, child_nodes, parent)
                else:
                    child_obj["case"] = "D"
                    child_obj["path"] = self.determine_path(child_obj, route, child_nodes)

        return child_nodes

    def determine_path(self, child_obj, route, child_nodes, parent=None):
        if not child_obj.get("case"):
            logging.warning("determine_path() - 'case' not in child_obj.")
            return False

        if child_obj["case"] == "A":
            return self.handle_case_a_path(route)
        elif child_obj["case"] == "B":
            return "1"
        elif child_obj["case"] == "C":
            return self.handle_case_c_path(route, child_nodes, parent)
        elif child_obj["case"] == "D":
            return self.handle_case_d_path(route)

    def handle_case_a_path(self, route):
        if route == "epoxide":
            return "1"
        elif route in self.cleaved_list:
            return "4" if "anhydride" in route else "5"
        elif "halogenated aliphatics" in route:
            return "2"
        return "3"

    def handle_case_c_path(self, route, child_nodes, parent):
        if route in self.op_esters:
            num_op_esters = self.count_op_esters(child_nodes)
            return "1" if num_op_esters > 4 else "2"
        elif "epoxide" in route:
            return "3"
        elif route in self.cleaved_list:
            return "5"
        elif parent is not None and any(mol in parent for mol in ["N", "P", "S", "O"]):
            logging.warning("SMILES contains N, S, P, or O. Returning qualitative value.")
            return "6"
        return "4"

    def handle_case_d_path(self, route):
        if "halogenated aliphatics" in route:
            return "1"
        elif "epoxide" in route:
            return "2"
        elif route in self.cleaved_list:
            return "3"
        return "4"

    def group_products(self, child_nodes):
        grouped_products = {}
        for child_obj in child_nodes:
            key = child_obj["case"] + child_obj["path"]
            grouped_products.setdefault(key, []).append(child_obj)
        return grouped_products

    def get_qsar_for_products_epi_api(self, parent, grouped_products):
        all_products_list = []

        for path_key, child_obj_list in grouped_products.items():
            case = path_key[0]
            path = path_key[1]
            route = child_obj_list[0]["routes"].lower()

            logging.info("Path key: {}\nRoute: {}\nUrl: {}".format(path_key, route, self.baseUrl))

            if path_key in ["A2", "C1", "D1"]:
                logging.info("Assigning qualitative values for case: {}, path: {}".format(case, path))
                child_obj_list = self.assign_qualitative_values(child_obj_list)
                all_products_list += child_obj_list
                continue

            response = requests.get(self.baseUrl, params={"smiles": parent})

            if response.status_code != 200:
                logging.warning("Error requesting half-life data. Status: {} Content: {}".format(
                    response.status_code, response.content))
                for child_obj in child_obj_list:
                    child_obj["error"] = "Error requesting half-life data from EPI."
                    child_obj["prop"] = "qsar"
                    child_obj["valid"] = False
                all_products_list += child_obj_list
                continue

            response_obj = json.loads(response.content)
            child_obj_list = self.handle_hl_response_epi_api(response_obj, parent, route, case, path, child_obj_list)
            all_products_list += child_obj_list

        return all_products_list

    def curate_api_response(self, response_obj, route: str = "") -> dict:
        """
        Parses EPI API response into the standard data shape.

        Now delegates to route-specific parsers ported from EPIReader.cs,
        replacing the previous mix of ad-hoc regex and halfLives iteration.

        The route parameter drives which C#-ported parser is used on the
        raw output text. Falls back to the halfLives JSON array if the
        text-based parse yields nothing (e.g. for routes the text parsers
        don't cover).
        """
        response_obj_new = {"data": []}

        hydrolysis_output_text = response_obj.get("hydrolysis", {}).get("output", "")

        logging.info("curate_api_response - route: '{}'\nOutput text:\n{}".format(
            route, hydrolysis_output_text))

        # --- Primary path: route-specific text parsing (ported from C#) ---
        if hydrolysis_output_text and route:
            parsed = parse_output_by_route(hydrolysis_output_text, route)
            if parsed:
                response_obj_new["data"] = parsed
                logging.info("curate_api_response - parsed {} result(s) via text parser.".format(len(parsed)))
                return response_obj_new

        # --- Fallback: halfLives JSON array (original Python logic) ---
        logging.info("curate_api_response - falling back to halfLives JSON.")
        hydrolysis_values = response_obj.get("hydrolysis", {})
        data_obj_template = {
            "chemical": None,
            "prop": None,
            "calc": "epi",
            "method": None,
            "data": None,
            "units": None,
        }

        for halflife_obj in hydrolysis_values.get("halfLives", []):
            if halflife_obj.get("ph") != 7.0:
                continue

            has_ka = halflife_obj.get("acidCatalyzed")
            has_kb = halflife_obj.get("baseCatalyzed")
            value = halflife_obj.get("value")
            units = halflife_obj.get("unit")

            if not value or value == 0:
                continue

            new_data_obj = dict(data_obj_template)
            new_data_obj["units"] = units
            new_data_obj["data"] = value

            if has_ka and not has_kb:
                new_data_obj["prop"] = "Ka"
                response_obj_new["data"].append(new_data_obj)
            elif has_kb and not has_ka:
                new_data_obj["prop"] = "Kb"
                response_obj_new["data"].append(new_data_obj)

        return response_obj_new

    def handle_hl_response_epi_api(self, response_obj, parent, route, case, path, child_obj_list):
        # Pass route through to curate_api_response so it can pick the right parser
        response_obj = self.curate_api_response(response_obj, route=route)
        logging.info("Curated response: {}".format(response_obj))
        return self._dispatch_hl_response(response_obj, parent, route, case, path, child_obj_list)

    def handle_hl_response(self, response_obj, parent, route, case, path, child_obj_list):
        return self._dispatch_hl_response(response_obj, parent, route, case, path, child_obj_list)

    def _dispatch_hl_response(self, response_obj, parent, route, case, path, child_obj_list):
        """
        Single dispatcher for both EPI API and legacy wrapper responses,
        replacing the duplicated if/elif blocks in the original two handle_hl_response methods.
        """
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
            return self.handle_case_d_results(path, response_obj, child_obj_list)
        elif case == "C":
            return self.handle_case_c_results(parent, path, route, response_obj, child_obj_list)
        return child_obj_list

    def handle_op_ester_values(self, response_obj, child_obj_list):
        op1_hl = None
        op2_hl = None

        for data_obj in response_obj.get("data", []):
            if data_obj["prop"] == "Kb":
                op1_hl = data_obj.get("data")
            elif data_obj["prop"] in ("Ka", "Kn"):
                op2_hl = data_obj.get("data")

        for child_obj in child_obj_list:
            route = child_obj.get("routes").lower()
            if route == self.op_esters[0]:
                child_obj["data"] = self.round_half_life(op1_hl)
            elif route == self.op_esters[1]:
                child_obj["data"] = self.round_half_life(op2_hl)

        return child_obj_list

    def handle_functional_group_case(self, route, parent, response_obj, child_obj_list):
        func_group = self.rdkit.get_functional_groups(route, parent)
        logging.warning("Functional groups: {}".format(func_group))
        if len(func_group) == 1:
            return self.hl_result_pattern("Kb", child_obj_list, response_obj)
        return self.assign_qualitative_values(child_obj_list)

    def handle_case_d_results(self, path, response_obj, child_obj_list):
        logging.info("Case: D{}".format(path))
        if path == "1":
            return self.assign_qualitative_values(child_obj_list)
        elif path == "2":
            return self.hl_result_pattern("Ka/n", child_obj_list, response_obj)
        elif path in ("3", "4"):
            return self.hl_result_pattern("Kb", child_obj_list, response_obj)

    def handle_case_c_results(self, parent, path, route, response_obj, child_obj_list):
        logging.info("Case: C{}".format(path))
        if path == "1":
            return self.assign_qualitative_values(child_obj_list)
        elif path == "2":
            return self.handle_op_ester_values(response_obj, child_obj_list)
        elif path == "3":
            return self.hl_result_pattern("Ka/n", child_obj_list, response_obj)
        elif path in ("4", "5"):
            return self.hl_result_pattern("Kb", child_obj_list, response_obj)
        elif path == "6":
            return self.assign_qualitative_values(child_obj_list)

    def hl_result_pattern(self, sort_prop, child_obj_list, response_obj):
        sorted_response = self.sort_k_by_atom_number(response_obj, sort_prop)
        if len(sorted_response) == 1:
            for child_obj in child_obj_list:
                child_obj["data"] = self.round_half_life(sorted_response[0]["data"])
            return child_obj_list
        elif len(sorted_response) > 1:
            return self.split_hl_values(child_obj_list, sorted_response, sort_prop)
        else:
            raise Exception("No HL values found for sort_prop: {}".format(sort_prop))

    def split_hl_values(self, child_obj_list, sorted_response, sort_prop="Kb"):
        num_products = len(child_obj_list)
        num_hls = len(sorted_response)
        mid_point = int(num_products / 2)

        logging.info("Splitting {} products across {} HLs.".format(num_products, num_hls))

        for child_obj in child_obj_list[0:mid_point]:
            child_obj["data"] = self.round_half_life(sorted_response[0]["data"])

        for child_obj in child_obj_list[mid_point:]:
            child_obj["data"] = self.round_half_life(sorted_response[1]["data"])

        return child_obj_list

    def make_qsar_request(self, request_dict):
        parent = request_dict.get("filtered_smiles")
        unique_schemes_count = int(request_dict.get("uniqueSchemesCount"))
        product_count = int(request_dict.get("productCount"))
        child_nodes = request_dict.get("childNodes")

        child_nodes = self.sort_products_by_case(parent, unique_schemes_count, product_count, child_nodes)
        grouped_products = self.group_products(child_nodes)

        logging.warning("Grouped Products: {}".format(grouped_products))

        return self.get_qsar_for_products_epi_api(parent, grouped_products)