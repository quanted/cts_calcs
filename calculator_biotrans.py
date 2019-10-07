import requests
import json
import time
import os
import logging
from bs4 import BeautifulSoup

from .calculator import Calculator
from .chemical_information import SMILESFilter



headers = {'Content-type': 'application/json', 'Accept': 'text/html'}



class BiotransCalc(Calculator):
	"""
	Biotransformer Calculator (http://biotransformer.ca)
	"""
	def __init__(self):
		Calculator.__init__(self)
		self.method = None
		self.name = "biotrans"
		self.baseUrl = os.environ.get('CTS_BIOTRANS_SERVER', "http://biotransformer.ca")
		self.query_url = self.baseUrl + '/queries.json'  # initiates request
		self.pred_url = self.baseUrl + '/queries/{}.json'  # returns status and results
		self.urlStruct = "/biotrans/rest/run"
		self.props = ["CYP450", "EC-BASED", "PHASEII", "HGUT", "ENVMICRO", "ALLHUMAN", "SUPERBIO"]
		self.biotrans_tasks = ["PREDICTION", "IDENTIFICATION"]
		self.request_timeout = 5
		self.max_response_wait = 30  # seconds (use request_timeout)
		self.gen_to_delay_map = {
			1: 1,
			2: 2,
			3: 5  # 3 gen limit, 5s delay
		}
		self.delay = 5
		self.data_received = False
		self.propMap = {}
		# POST object for biotransformer API:
		self.postData = {
			"biotransformer_option": self.props[0],
			"number_of_steps": 1,
			"query_input": "",  # the smiles
			"task_type": self.biotrans_tasks[0]
		}
		# POST object for CTS API biotransformer calc:
		self.api_post = {
			"chemical": "",
			"prop": "",
			"gen_limit": 1
		}
		self.meta_info = {
			'metaInfo': {
				'model': "biotrans",
				'collection': "qed",
				'modelVersion': "",
				'description': "BioTransformer is a software tool that predicts small molecule metabolism in mammals, their gut microbiota, as well as the soil/aquatic microbiota.",
				'timestamp': self.gen_jid(),
				'url': "http://biotransformer.ca/",
				'props': self.props,
			}
		}

	def current_milli_time(self):
		"""
		Returns time in milliseconds.
		"""
		return int(round(time.time() * 1000))

	def get_query_id(self, api_response):
		"""
		Creates object from HTML, then gets query id from the 'data-query-id' div attribute.
		"""
		soup = BeautifulSoup(api_response.content, features="html.parser")  # creates an object from the html response
		content = soup.find("div", {"id": "query-status"})  # gets div from response that has query/response id
		request_id = content.get('data-query-id')  # gets 'data-query-id' attr from div, which is the request's id
		return request_id

	def curate_data(self, metabolites_list):
		"""
		Creates a nested dictionary from the biotransformer predictions data.
		"""
		has_parent = set()
		all_items = {}
		met_id = 1
		for metabolite_obj in metabolites_list:
			parent_obj = metabolite_obj['substrates'][0]
			child_obj = metabolite_obj['products'][0]
			parent = parent_obj['smiles']
			child = child_obj['smiles']
			if parent not in all_items:
				all_items[parent] = {}
				all_items[parent]['id'] = met_id
				all_items[parent]['name'] = "<img class='blank_node' src='/static_qed/cts/images/loader_node.gif' />"
				# all_items[parent]['html'] = parent_obj
				all_items[parent]['data'] = parent_obj
				all_items[parent]['children'] = []
				met_id += 1
			if child not in all_items:
				all_items[child] = {}
				all_items[child]['id'] = met_id
				all_items[child]['name'] = "<img class='blank_node' src='/static_qed/cts/images/loader_node.gif' />"
				# all_items[child]['html'] = child_obj
				all_items[child]['data'] = child_obj
				all_items[child]['children'] = []
				met_id += 1
			all_items[parent]['children'].append(all_items[child])
			has_parent.add(child)
		result = {}
		for key, value in all_items.items():
			if key not in has_parent:
				# result[key] = value
				result['tree'] = value
		return result

	def make_query_request(self, api_query):
		"""
		Makes request to initiate prediction calculations.
		Returns an ID for looking up response and status.
		"""
		try:
			return requests.post(self.query_url, json=api_query, headers=headers, timeout=self.request_timeout)
		except Exception as e:
			logging.warning("Exception in calculator_biotrans: {}".format(e))
			return None

	def make_status_request(self, query_id):
		"""
		Makes request for predictions.
		"""
		try:
			result = requests.get(self.pred_url.format(query_id))  # /queries/[id].json
			return json.loads(result.content)
		except Exception as e:
			logging.warning("Exception in calculator_biotrans: {}".format(e))
			return None

	def poll_predictions(self, query_id, gen_limit):
		"""
		Polls query id endpoint for predictions data.
		"""
		delay_sum = 0
		delay = self.gen_to_delay_map[gen_limit]
		while True:
			logging.info("Sleeping {}s before making request".format(delay))
			time.sleep(delay)
			result_obj = self.make_status_request(query_id)
			status = result_obj.get('status')
			logging.info("Query Status: {}".format(status))
			if status == "Done":
				return result_obj
			elif status == "failed":
				logging.warning("calculator_biotrans - received 'failed' from API")
				raise Exception("Cannot reach biotransformer")
			elif float(delay_sum) / float(self.max_response_wait) > 1.0:
				logging.warning("Exceeded max wait for biotransformer response ({}s)".format(self.max_response_wait))
				raise Exception("Exceeded max tries")
			delay_sum += delay

	def data_request_handler(self, request_dict):
		"""
		Entrypoint for handling biotransformer data request.
		"""
		chemical = request_dict['chemical']
		# prop = request_dict['prop']
		prop = request_dict.get('prop', self.props[0])
		gen_limit = int(request_dict.get('gen_limit', 1))

		if 'metabolizer_post' in request_dict:
			prop = request_dict['metabolizer_post'].get('prop', self.props[0])
			gen_limit = int(request_dict['metabolizer_post'].get('gen_limit', 1))

		if not prop in self.props:
			return {'status': False, 'error': "Select an available prop: {}".format(self.props)}

		start_time = self.current_milli_time()
		
		api_query = dict(self.postData)
		api_query['biotransformer_option'] = prop
		api_query['query_input'] = chemical
		api_query['number_of_steps'] = gen_limit

		api_response = self.make_query_request(api_query)
		if not api_response:
			# error occurred requesting data
			return  # TODO: more granular error handling

		query_id = self.get_query_id(api_response)
		response = self.poll_predictions(query_id, gen_limit)

		if not response or not 'predictions' in response:
			raise Exception({'error': "Could not retrieve predictions."})

		metabolites_list = response['predictions'][0]['biotransformations']
		tree_obj = self.curate_data(metabolites_list)

		tree_obj['total_products'] = response['predictions'][0]['nr_of_biotransformations']

		if tree_obj['total_products'] < 1:
			raise Exception("No products.")

		execution_time = self.current_milli_time() - start_time
		logging.info("Execution time (ms): {}".format(execution_time))

		_response_obj = {
            'calc': "chemaxon",  # todo: change to metabolizer, change in template too
            'prop': "products",
            'node': request_dict.get('node'),
            'data': tree_obj['tree'],
            'total_products': tree_obj['total_products'],
            'chemical': request_dict.get('chemical'),
            'workflow': 'gentrans',
            'run_type': 'batch',
            'request_post': request_dict       
        }

		return _response_obj