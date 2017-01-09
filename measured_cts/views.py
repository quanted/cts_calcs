import logging
import os
import requests
import json
from measured_calculator import MeasuredCalc
import redis
from cts_app.cts_calcs.smilesfilter import parseSmilesByCalculator


def request_manager(request):

	calc = request.POST.get("calc")
	# props = request.POST.getlist("props[]")
	try:
		props = request.POST.get("props[]")
		if not props:
			props = request.POST.getlist("props")
	except AttributeError:
		props = request.POST.get("props")

	structure = request.POST.get("chemical")
	sessionid = request.POST.get('sessionid')
	node = request.POST.get('node')
	run_type = request.POST.get('run_type')
	prop = request.POST.get('prop')

	logging.info("Incoming data to Measured: {}, {}, {} (calc, props, chemical)".format(calc, props, structure))

	post_data = {
		"calc": calc,
		"props": props,
	    'node': node
	}

	logging.info("Measured receiving SMILES: {}".format(structure))

	try:
		filtered_smiles = parseSmilesByCalculator(structure, "epi") # call smilesfilter
		if '[' in filtered_smiles or ']' in filtered_smiles:
			logging.warning("EPI ignoring request due to brackets in SMILES..")
			post_data.update({'error': "EPI Suite cannot process charged species or metals (e.g., [S+], [c+])"})
			return HttpResponse(json.dumps(post_data), content_type='application/json')
	except Exception as err:
		logging.warning("Error filtering SMILES: {}".format(err))
		post_data.update({'error': "Cannot filter SMILES for Measured data"})
		return HttpResponse(json.dumps(post_data), content_type='application/json')

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
			data_obj.update({'node': node, 'request_post': request.POST})

			# push one result at a time if node/redis:
			result_json = json.dumps(data_obj)
			return HttpResponse(result_json, content_type='application/json')

	except Exception as err:
		logging.warning("Exception occurred getting Measured data: {}".format(err))
		for prop in props:
			data_obj = {
				'data': "error - cannot find measured data for structure",
				'calc': "measured",
				'prop': prop,
				'node': node,
				'request_post': request.POST
			}
			return HttpResponse(json.dumps(data_obj), content_type='application/json')