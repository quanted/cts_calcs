import logging
import os
import requests
import json
from test_calculator import TestCalc
try:
	from cts_app.cts_calcs.smilesfilter import parseSmilesByCalculator
except ImportError as e:
	from cts_calcs.smilesfilter import parseSmilesByCalculator

import redis


redis_hostname = os.environ.get('REDIS_HOSTNAME')
redis_conn = redis.StrictRedis(host=redis_hostname, port=6379, db=0)


def request_manager(request):
	"""
	less_simple_proxy takes a request and
	makes the proper call to the TEST web services.
	it relies on the test_calculator to do such.

	input: {"calc": [calculator], "prop": [property]}
	output: returns data from TEST server
	"""

	TEST_URL = os.environ["CTS_TEST_SERVER"]
	postData = {}

	calc = request.POST.get("calc")
	# prop = request.POST.get("prop")

	try:
		props = request.POST.get("props[]")
		if not props:
			props = request.POST.getlist("props")
	except AttributeError:
		props = request.POST.get("props")

	calc_data = request.POST.get('checkedCalcsAndProps')
	chemical = request.POST.get("chemical")
	sessionid = request.POST.get('sessionid')
	node = request.POST.get('node')
	mass = request.POST.get('mass')  # for water solubility
	run_type = request.POST.get('run_type')
    workflow = request.POST.get('workflow')

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
		redis_conn.publish(sessionid, json.dumps(postData))
	# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	logging.info("TEST Filtered SMILES: {}".format(filtered_smiles))

	calcObj = TestCalc()
	test_results = []
	for prop in props:

		data_obj = {
			'chemical': chemical,
			'calc': calc,
			'prop':prop,
			'node': node,
			'request_post': request.POST,
			'run_type': run_type,
			'workflow': workflow,
		}

		try:
			logging.info("Calling TEST for {} data...".format(prop))

			response = calcObj.makeDataRequest(filtered_smiles, calc, prop)
			response_json = json.loads(response.content)

			logging.info("TEST response data for {}: {}".format(prop, response_json))

			# sometimes TEST data successfully returns but with an error:
			if response.status_code != 200:
				postData['data'] = "TEST could not process chemical"
			else:
				test_data = response_json['properties'][calcObj.propMap[prop]['urlKey']]
				if test_data == -9999:
					data_obj['data'] = "N/A"
				elif prop == 'water_sol':
					data_obj['data'] = calcObj.convertWaterSolubility(mass, test_data)
				else:
					data_obj['data'] = test_data
				
			result_json = json.dumps(data_obj)
			redis_conn.publish(sessionid, result_json)

		except Exception as err:
			logging.warning("Exception occurred getting TEST data: {}".format(err))
			data_obj.update({'data': "timed out", 'request_post': request.POST})

			logging.info("##### session id: {}".format(sessionid))

			redis_conn.publish(sessionid, json.dumps(data_obj))