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
		self.db = None  # opens cts database (set in connection function)
		self.chem_info_collection = None  # chem info data collection (set in connection function)
		# self.pchem_collection = self.db.pchem  # pchem data collection
		# self.products_collection = self.db.products  # transformation products collection
		self.db_conn_timeout = 1
		self.is_connected = False

		# Keys for chem info collection document entry:
		self.chem_info_keys = ["chemical", "orig_smiles", "smiles",
			"preferredName", "iupac", "formula", "casrn", "cas",
			"dsstoxSubstanceId", "mass", "exactmass", "structureData"]

		# Keys for pchem collection document entry:
		# self.pchem_keys = ["smiles", "calc", "prop", "data", "method"]

	def connect_to_db(self):
		"""
		Tries to connect to mongodb.
		"""
		try:
			self.mongodb_conn = pymongo.MongoClient(port=27017, serverSelectionTimeoutMS=1000, connectTimeoutMS=1000)
			self.is_connected = True
			self.db = self.mongodb_conn.cts  # opens cts database
			self.chem_info_collection = self.db.chem_info  # chem info data collection
			self.test_db_connection()
		except pymongo.errors.ConnectionFailure as e:
			logging.warning("Unable to connect to db: {}".format(e))
			self.is_connected = False
		except pymongo.errors.ConfigurationError as e:
			logging.warning("Config error setting up mongodb client: {}".format(e))
			logging.warning("Unable to connect to db.")
			self.is_connected = False

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
		new_query_obj = dict()
		for key, val in query_obj.items():
			if key in self.chem_info_keys:
				new_query_obj[key] = val
		return new_query_obj

	# def find_cheminfo_with_dsstox(self, user_chemical):
	# 	"""
	# 	Gets DSSTOX to use as key for chem-info DB objects.
	# 	"""
	# 	if not self.is_connected:
	# 		return None


	def find_chem_info_document(self, query_obj):
		"""
		Searches chem info collection for document matching chemical.
		Returns chem info data if it exists, or None if it doesn't.
		"""
		if not self.is_connected:
			return None
		chem_info_result = self.chem_info_collection.find_one(query_obj)  # searches db
		return chem_info_result

	def insert_chem_info_data(self, molecule_obj):
		"""
		Inserts chem info data into chem info collection.
		Returns document unique _id.
		"""
		if not self.is_connected and not molecule_obj:
			return None
		db_object = self.create_query_obj(molecule_obj)
		chem_info_obj = self.chem_info_collection.insert_one(db_object)  # inserts query object
		return chem_info_obj

	# def insert_pchem_data(self, pchem_data):
	# 	"""
	# 	"""