import requests
import json
import logging
import os
try:
    from cts_app.cts_calcs.calculator import Calculator
except ImportError as e:
    from cts_calcs.calculator import Calculator

headers = {'Content-Type': 'application/json'}


class TestCalc(Calculator):
    """
	TEST Suite Calculator
	"""

    def __init__(self):

        Calculator.__init__(self)

        self.postData = {"smiles" : ""}
        self.name = "test"
        self.baseUrl = os.environ['CTS_TEST_SERVER']
        self.urlStruct = "/api/TEST/{}/{}" # method, property

        # self.methods = ['hierarchical']
        self.methods = ['FDAMethod', 'HierarchicalMethod']
        # map workflow parameters to test
        self.propMap = {
            'melting_point': {
               'urlKey': 'MeltingPoint'
            },
            'boiling_point': {
               'urlKey': 'BoilingPoint'
            },
            'water_sol': {
               'urlKey': 'WaterSolubility'
            },
            'vapor_press': {
               'urlKey': 'VaporPressure'
            }
            # 'henrys_law_con': ,
            # 'kow_no_ph': 
        }


    def getPostData(self, calc, prop, method=None):
        return {"identifiers": {"SMILES": ""}}


    def makeDataRequest(self, structure, calc, prop, method=None):
        post = self.getPostData(calc, prop)
        post['identifiers']['SMILES'] = structure # set smiles
        test_prop = self.propMap[prop]['urlKey'] # prop name TEST understands
        url = self.baseUrl + self.urlStruct.format('FDAMethod', test_prop)
        try:
            response = requests.post(url, data=json.dumps(post), headers=headers, timeout=60)
        except requests.exceptions.ConnectionError as ce:
            logging.info("connection exception: {}".format(ce))
            # return None
            raise
        except requests.exceptions.Timeout as te:
            logging.info("timeout exception: {}".format(te))
            # return None
            raise

        self.results = response
        return response


    def convertWaterSolubility(self, mass, test_datum):
        """
        Converts water solubility from log(mol/L) => mg/L
        """
        return 1000 * float(mass) * 10**-(test_datum)


    def data_request_handler(self, request_dict):
        
        TEST_URL = os.environ["CTS_TEST_SERVER"]
        postData = {}

        calc = request.get("calc")
        # prop = request.get("prop")

        try:
            props = request.get("props[]")
            if not props:
                props = request.getlist("props")
        except AttributeError:
            props = request.get("props")

        calc_data = request.get('checkedCalcsAndProps')
        chemical = request.get("chemical")
        sessionid = request.get('sessionid')
        node = request.get('node')
        mass = request.get('mass')  # for water solubility
        run_type = request.get('run_type')
        workflow = request.get('workflow')

        if calc_data:
            calc = "test"
            props = calc_data['test']  # list of props

        postData = {
            "chemical": chemical,
            "calc": calc,
            'run_type': run_type,
            'workflow': workflow,
            # "prop": prop
            # "props": props
        }

        # filter smiles before sending to TEST:
        # ++++++++++++++++++++++++ smiles filtering!!! ++++++++++++++++++++
        try:
            filtered_smiles = parseSmilesByCalculator(chemical, calc) # call smilesfilter
        except Exception as err:
            logging.warning("Error filtering SMILES: {}".format(err))
            postData.update({'data': "Cannot filter SMILES for TEST data"})
            self.redis_conn.publish(sessionid, json.dumps(postData))
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        logging.info("TEST Filtered SMILES: {}".format(filtered_smiles))

        self = TestCalc()
        test_results = []
        for prop in props:

            data_obj = {
                'chemical': filtered_smiles,
                'calc': calc,
                'prop':prop,
                'node': node,
                'request_post': request.POST,
                'run_type': run_type,
                'workflow': workflow,
            }

            try:
                logging.info("Calling TEST for {} data...".format(prop))

                response = self.makeDataRequest(filtered_smiles, calc, prop)
                response_json = json.loads(response.content)

                logging.info("TEST response data for {}: {}".format(prop, response_json))

                # sometimes TEST data successfully returns but with an error:
                if response.status_code != 200:
                    postData['data'] = "TEST could not process chemical"
                else:
                    test_data = response_json['properties'][self.propMap[prop]['urlKey']]
                    if test_data == -9999:
                        data_obj['data'] = "N/A"
                    elif prop == 'water_sol':
                        data_obj['data'] = self.convertWaterSolubility(mass, test_data)
                    else:
                        data_obj['data'] = test_data
                    
                result_json = json.dumps(data_obj)
                self.redis_conn.publish(sessionid, result_json)

            except Exception as err:
                logging.warning("Exception occurred getting TEST data: {}".format(err))
                data_obj.update({'data': "timed out", 'request_post': request.POST})

                logging.info("##### session id: {}".format(sessionid))

                self.redis_conn.publish(sessionid, json.dumps(data_obj))