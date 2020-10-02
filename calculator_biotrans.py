import requests
import json
import time
import os
import logging
from bs4 import BeautifulSoup

from .calculator import Calculator
from .chemical_information import SMILESFilter



headers = {"Content-type": "application/json", "Accept": "text/html"}



class BiotransCalc(Calculator):
	"""
	Biotransformer Calculator (local wrapper).
	"""
	def __init__(self):
		Calculator.__init__(self)
		self.method = None
		self.name = "biotrans"
		self.baseUrl = os.environ.get("CTS_BIOTRANS_SERVER")
		self.urlStruct = "/bt/rest/run"
		self.biotrans_api_url = self.baseUrl + self.urlStruct
		self.props = ["ecbased", "cyp450", "phaseII", "hgut", "superbio", "envimicro"]
		self.biotrans_tasks = ["pred", "cid"]
		self.request_timeout = 30
		self.propMap = {}

		self.meta_info = {
			"metaInfo": {
				"model": "biotrans",
				"collection": "qed",
				"modelVersion": "",
				"description": "BioTransformer is a software tool that predicts small molecule metabolism in mammals, their gut microbiota, as well as the soil/aquatic microbiota.",
				"timestamp": self.gen_jid(),
				"url": "http://biotransformer.ca/",
				"props": self.props,
			}
		}

	def make_request(self, url, post_data):
		"""
		Makes request to initiate prediction calculations.
		Returns an ID for looking up response and status.
		"""
		try:
			return requests.post(url, json=post_data, headers=headers, timeout=self.request_timeout)
		except Exception as e:
			logging.warning("Exception in calculator_biotrans: {}".format(e))
			return {"error": "Error making request to biotransformer."}

	def data_request_handler(self, request_dict):
		"""
		Entrypoint for handling biotransformer data request.
		"""

		metabolizer_data = request_dict.get("metabolizer_post")
		chemical = request_dict["chemical"]
		prop = metabolizer_data.get("prop")
		gen_limit = int(request_dict.get("gen_limit", 1))

		if not prop in self.props:
			return {"status": False, "error": "Select an available prop: {}".format(self.props)}

		post_data = {
			"smiles": chemical,
			"prop": prop,
			"gen_limit": gen_limit
		}

		response = self.make_request(self.biotrans_api_url, post_data)  # return tree and total_products

		if "error" in response:
			return {"status": False, "error": response["error"]}

		response_obj = json.loads(response.content.decode("utf-8"))

		return {
            "calc": "biotrans",  # todo: change to metabolizer, change in template too
            "prop": "products",
            "node": request_dict.get("node"),
            "data": response_obj["data"]["tree"],
            "total_products": response_obj["data"]["total_products"],
            "chemical": request_dict.get("chemical"),
            "workflow": request_dict.get("workflow"),
			"run_type": request_dict.get("run_type"),
            "request_post": request_dict       
        }
