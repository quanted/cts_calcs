"""
Handles CTS mongodb interactions.
"""

import pymongo
import datetime
import pytz
import logging
import os



class MongoDBHandler:

	def __init__(self):
		# MongoDB Settings:
		self.db = None  # opens cts database (set in connection function)
		self.chem_info_collection = None  # chem info data collection (set in connection function)
		self.pchem_collection = None  # pchem data collection
		self.db_conn_timeout = 1
		self.is_connected = False
		self.mongodb_conn = None
		self.mongodb_host = os.environ.get('CTS_DB_HOST')

		# Keys for chem info collection document entry:
		self.chem_info_keys = ["dsstoxSubstanceId", "chemical", "orig_smiles",
			"smiles", "preferredName", "iupac", "formula", "casrn",
			"cas", "mass", "exactMass", "structureData"]
		self.extra_chem_info_Keys = ["node_image", "popup_image"]  # html wrappers w/ images for product nodes and popups

		# Keys for pchem collection document entry:
		self.pchem_keys = ["dsstoxSubstanceId", "calc", "prop", "data", "method", "ph"]

		# Keys for dtxcid -> dtxsid document entries (see DSST_IDs.csv):
		self.dtxcid_keys = ["DTXCID", "DTXSID", "CASRN", "PreferredName"]

	def connect_to_db(self):
		"""
		Tries to connect to mongodb.
		"""
		try:
			logging.info("(mongodb_handler.py) Connecting to MongoDB at: {}".format(self.mongodb_host	))
			self.mongodb_conn = pymongo.MongoClient(host=self.mongodb_host, serverSelectionTimeoutMS=200, connectTimeoutMS=200)
			self.is_connected = True
			self.db = self.mongodb_conn.cts  # opens cts database
			# self.chem_info_collection = self.db.chem_info  # chem info data collection
			self.pchem_collection = self.db.pchem  # pchem data collection
			self.dtxcid_collection = self.db.dtxcid  # dtxcid data collection
			self.test_db_connection()
		except pymongo.errors.ConnectionFailure as e:
			logging.warning("(mongodb_handler.py) Unable to connect to db: {}".format(e))
			self.is_connected = False
			self.mongodb_conn.close()
		except pymongo.errors.ConfigurationError as e:
			logging.warning("Config error setting up mongodb client: {}".format(e))
			logging.warning("(mongodb_handler.py) Unable to connect to db.")
			self.is_connected = False
			self.mongodb_conn.close()
		except Exception as e:
			logging.warning("(mongodb_handler.py) Error connecting to db: {}".format(e))
			self.is_connected = False
			self.mongodb_conn.close()

	def test_db_connection(self):
		"""
		Tests mongodb connection, sets is_connected to True
		if connected.
		"""
		try:
			self.mongodb_conn.server_info()
			self.is_connected = True
		except pymongo.errors.ServerSelectionTimeoutError as e:
			logging.warning("(mongodb_handler.py) Unable to connect to db.")
			self.is_connected = False

	def gen_jid(self):
		ts = datetime.datetime.now(pytz.UTC)
		localDatetime = ts.astimezone(pytz.timezone('US/Eastern'))
		jid = localDatetime.strftime('%Y%m%d%H%M%S%f')
		return jid

	# def create_chem_info_document(self, query_obj):
	# 	"""
	# 	Creates chem info object for querying and inserting.
	# 	"""
	# 	new_query_obj = dict()
	# 	for key, val in query_obj.items():
	# 		if key in self.chem_info_keys + self.extra_chem_info_Keys:
	# 			new_query_obj[key] = val
	# 	return new_query_obj

	# def find_chem_info_document(self, query_obj):
	# 	"""
	# 	Searches chem info collection for document matching chemical.
	# 	Returns chem info data if it exists, or None if it doesn't.
	# 	"""
	# 	if not self.is_connected:
	# 		return None
	# 	chem_info_result = self.chem_info_collection.find_one(query_obj)  # searches db
	# 	return chem_info_result

	# def insert_chem_info_data(self, molecule_obj):
	# 	"""
	# 	Inserts chem info data into chem info collection.
	# 	Returns document unique _id.
	# 	"""
	# 	if not self.is_connected or not molecule_obj:
	# 		return None
	# 	db_object = self.create_chem_info_document(molecule_obj)
	# 	chem_info_obj = self.chem_info_collection.insert_one(db_object)  # inserts query object
	# 	return chem_info_obj

	def create_pchem_document(self, query_obj):
		"""
		Creates p-chem object for querying and inserting.
		"""
		if not self.is_connected or not query_obj:
			return None
		new_query_obj = dict()
		for key, val in query_obj.items():
			if key in self.pchem_keys:
				new_query_obj[key] = val
		return new_query_obj

	def find_pchem_document(self, query_obj):
		"""
		Searches pchem collection for document matching chemical.
		Returns pchem data if it exists, or None if it doesn't.
		"""
		# if not self.is_connected or not query_obj:
		# 	return None
		pchem_result = self.pchem_collection.find_one(query_obj)  # searches db
		return pchem_result

	# def insert_pchem_data(self, pchem_obj):
	# 	"""
	# 	Inserts pchem data into chem info collection.
	# 	Returns document unique _id.
	# 	"""
	# 	# if not self.is_connected or not pchem_obj:
	# 	# 	return None
	# 	print("Inserting {} into db.".format(pchem_obj))
	# 	db_object = self.create_pchem_document(pchem_obj)
	# 	pchem_obj = self.pchem_collection.insert_one(db_object)  # inserts query object
	# 	return pchem_obj

	# def create_dtxcid_document(self, query_obj):
	# 	"""
	# 	Creates dtxcid object for querying and inserting.
	# 	"""
	# 	if not self.is_connected or not query_obj:
	# 		return None
	# 	new_query_obj = dict()
	# 	for key, val in query_obj.items():
	# 		if key in self.dtxcid_keys:
	# 			new_query_obj[key] = val
	# 	return new_query_obj

	def find_dtxcid_document(self, query_obj):
		"""
		Searches dtxcid collection for document matching chemical.
		Returns dtxcid data if it exists, or None if it doesn't.
		"""
		# if not self.is_connected or not query_obj:
		# 	return None
		dtxcid_result = self.dtxcid_collection.find_one(query_obj)  # searches db
		return dtxcid_result

	# def insert_dtxcid_data(self, dtxcid_obj):
	# 	"""
	# 	Inserts dtxcid data into dtxcid collection.
	# 	Returns document unique _id.
	# 	"""
	# 	# if not self.is_connected or not dtxcid_obj:
	# 	# 	return None
	# 	print("Inserting {} into db.".format(dtxcid_obj))
	# 	db_object = self.create_dtxcid_document(dtxcid_obj)
	# 	dtxcid_obj = self.dtxcid_collection.insert_one(db_object)  # inserts query object
	# 	return dtxcid_obj