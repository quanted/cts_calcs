import logging
import os
import requests
import json
from measured_calculator import MeasuredCalc
import redis
try:
	from cts_app.cts_calcs.smilesfilter import parseSmilesByCalculator
except ImportError as e:
	from cts_calcs.smilesfilter import parseSmilesByCalculator


redis_hostname = os.environ.get('REDIS_HOSTNAME')
redis_conn = redis.StrictRedis(host=redis_hostname, port=6379, db=0)


def request_manager(request):

	calc = request.POST.get("calc")
	# props = request.POST.getlist("props[]")
	try:
		props = request.POST.get("props[]")
		if not props:
			props = request.POST.getlist("props")
	except AttributeError:
		props = request.POST.get("props")

	chemical = request.POST.get("chemical")
	sessionid = request.POST.get('sessionid')
	node = request.POST.get('node')
	run_type = request.POST.get('run_type')
	prop = request.POST.get('prop')

	logging.info("Incoming data to Measured: {}, {}, {} (calc, props, chemical)".format(calc, props, chemical))

	post_data = {
		"chemical": chemical,
		"calc": calc,
		"props": props,
	    'node': node,
	    'run_type': run_type,
		'workflow': workflow,
	}

	logging.info("Measured receiving SMILES: {}".format(chemical))

	try:
		filtered_smiles = parseSmilesByCalculator(chemical, "epi") # call smilesfilter
		if '[' in filtered_smiles or ']' in filtered_smiles:
			logging.warning("EPI ignoring request due to brackets in SMILES..")
			post_data.update({'error': "EPI Suite cannot process charged species or metals (e.g., [S+], [c+])"})
			redis_conn.publish(sessionid, json.dumps(post_data))
	except Exception as err:
		logging.warning("Error filtering SMILES: {}".format(err))
		post_data.update({'error': "Cannot filter SMILES for Measured data"})
		redis_conn.publish(sessionid, json.dumps(post_data))

	logging.info("Measured Filtered SMILES: {}".format(filtered_smiles))

	calcObj = MeasuredCalc()

	try:
		response = calcObj.makeDataRequest(filtered_smiles) # make call for data!
		measured_data = json.loads(response.content)

		if run_type == 'rest':
			props = [prop]
		
		# get requested properties from results:
		for prop in props:
			data_obj = calcObj.getPropertyValue(prop, measured_data)
			data_obj.update({
				'node': node, 
				'request_post': request.POST,
				'chemical': chemical,
				'run_type': run_type,
				'workflow': workflow,
			})

			# push one result at a time if node/redis:
			result_json = json.dumps(data_obj)
			redis_conn.publish(sessionid, result_json)

	except Exception as err:
		logging.warning("Exception occurred getting Measured data: {}".format(err))
		for prop in props:
			data_obj = {
				'data': "error - cannot find measured data for chemical",
				'calc': "measured",
				'prop': prop,
				'node': node,
				'request_post': request.POST,
				'chemical': chemical,
				'run_type': run_type,
				'workflow': workflow,
			}

			redis_conn.publish(sessionid, json.dumps(data_obj))