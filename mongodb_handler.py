"""
Handles CTS mongodb interactions.
"""

import pymongo
import datetime
import pytz
import logging



class MongoDBHandler:

	def __init__(self):
		# MongoDB Settings:
		self.mongodb_conn = pymongo.MongoClient(port=27017)
		self.db = self.mongodb_conn.cts  # opens cts database
		self.chem_info_collection = self.db.chem_info  # chem info data collection
		# self.pchem_collection = self.db.pchem  # pchem data collection
		# self.products_collection = self.db.products  # transformation products collection
		self.db_conn_timeout = 1

		# Keys for chem info collection document entry:
		self.chem_info_keys = ["chemical", "orig_smiles", "smiles", "preferredName", 
			"iupac", "formula", "casrn", "cas", "dsstoxSubstanceId", "mass", "exactmass"]

		# Keys for pchem collection document entry:
		# self.pchem_keys = ["smiles", "calc", "prop", "data", "method"]

		self.test_db_connection()

	def test_db_connection(self):
		"""
		Tests mongodb connection, sets is_connected to True
		if connected.
		"""
		try:
			self.mongodb_conn.server_info()
			self.is_connected = True
		except pymongo.errors.ServerSelectionTimeoutError as e:
			logging.warning("Unable to connect to db.")
			self.is_connected = False

	def gen_jid(self):
		ts = datetime.datetime.now(pytz.UTC)
		localDatetime = ts.astimezone(pytz.timezone('US/Eastern'))
		jid = localDatetime.strftime('%Y%m%d%H%M%S%f')
		return jid

	def create_query_obj(self, query_obj):
		"""
		Creates chem info object for querying and inserting.
		"""
		query_obj = dict()
		for key, val in query_obj.items():
			if key in self.chem_info_keys:
				query_obj[key] = val
		return query_obj

	def find_chem_info_document(self, query_obj, is_connected=False):
		"""
		Searches chem info collection for document matching chemical.
		Returns chem info data if it exists, or None if it doesn't.
		"""
		if not is_connected:
			return None
		chem_info_result = self.chem_info_collection.find_one(query_obj)  # searches db
		return chem_info_result

	def insert_chem_info_data(self, molecule_obj, is_connected=False):
		"""
		Inserts chem info data into chem info collection.
		Returns document unique _id.
		"""
		if not is_connected:
			return None
		chem_info_obj = self.chem_info_collection.insert_one(molecule_obj)  # inserts query object
		return chem_info_obj

	# def insert_pchem_data(self, pchem_data):
	# 	"""
	# 	"""