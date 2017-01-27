import logging
import requests
import json
import redis
import os

from sparc_calculator import SparcCalc

from cts_calcs.epi_cts import views as epi_views
from cts_calcs.measured_cts import views as measured_views


redis_hostname = os.environ.get('REDIS_HOSTNAME')
redis_conn = redis.StrictRedis(host=redis_hostname, port=6379, db=0)


def request_manager(request):

    calc = request.POST.get("calc")
    try:
        props = request.POST.get("props[]")
        if not props:
            props = request.POST.getlist("props")
    except AttributeError:
        props = request.POST.get("props")

    structure = request.POST.get("chemical")
    sessionid = request.POST.get('sessionid')
    node = request.POST.get('node')
    prop = request.POST.get('prop')
    run_type = request.POST.get('run_type')

    # logging.info("Incoming data to SPARC: {}, {}, {} (calc, props, chemical)".format(calc, props, structure))

    post_data = {
        "calc": "sparc",
        "props": props,
        'node': node
    }

    ############### Filtered SMILES stuff!!! ###########################
    # filtered_smiles = parseSmilesByCalculator(structure, "sparc") # call smilesfilter
    ###########################################################################################




    # # need to get melting point for sparc calculations.
    # # Try Measured, but what about EPI? I wouldn't mess with TEST..
    # melting_point_request = {
    #     'calc': "measured",  # should prob be measured
    #     'props': ['melting_point'],
    #     'chemical': structure,
    #     'sessionid': sessionid
    # }
    # request = NotDjangoRequest(melting_point_request)
    # melting_point_response = measured_views.request_manager(request)
    # logging.warning("MELTING POINT RESPONSE: {}".format(melting_point_response))
    # logging.warning("MELTING POINT RESPONSE TYPE: {}".format(type(melting_point_response)))




    calcObj = SparcCalc(structure)

    ion_con_response, kow_wph_response, multi_response = None, None, None
    sparc_results = []

    # synchronous calls for ion_con, kow_wph, and the rest:

    # don't need this loop, just do "if 'ion_con' in prop: make request"

    logging.warning("sparc props: {}".format(props))

    sparc_response = {
        'calc': 'sparc',
        # 'prop': prop,
        'node': node,
        'chemical': structure,
        'request_post': request.POST
    }

    if run_type == 'rest':
        props = [prop]


    try:
        if 'ion_con' in props:
            response = calcObj.makeCallForPka() # response as d ict returned..
            pka_data = calcObj.getPkaResults(response)
            sparc_response.update({'data': pka_data, 'prop': 'ion_con'})
            logging.info("response: {}".format(sparc_response))
            result_json = json.dumps(sparc_response)
            redis_conn.publish(sessionid, result_json)

        if 'kow_wph' in props:
            ph = request.POST.get('ph') # get PH for logD calculation..
            response = calcObj.makeCallForLogD() # response as dict returned..
            sparc_response.update({'data': calcObj.getLogDForPH(response, ph), 'prop': 'kow_wph'})
            logging.info("response: {}".format(sparc_response))
            result_json = json.dumps(sparc_response)
            redis_conn.publish(sessionid, result_json)

        multi_response = calcObj.makeDataRequest()
        if 'calculationResults' in multi_response:
            multi_response = calcObj.parseMultiPropResponse(multi_response['calculationResults'])
            for prop_obj in multi_response:
                if prop_obj['prop'] in props and prop_obj != 'ion_con' and prop_obj != 'kow_wph':
                    prop_obj.update({'node': node, 'chemical': structure, 'request_post': request.POST})
                    logging.info("response: {}".format(prop_obj))
                    result_json = json.dumps(prop_obj) 
                    redis_conn.publish(sessionid, result_json)

    except Exception as err:
        logging.warning("Exception occurred getting SPARC data: {}".format(err))

        for prop in props:

            post_data.update({
                'data': "data request timed out",
                'prop': prop,
                'request_post': request.POST
            })

            redis_conn.publish(sessionid, json.dumps(post_data))

# patch for freeing celery from django while calc views
# are still relying on django.http Request...
class NotDjangoRequest(object):
    def __init__(self, post_obj):
        self.POST = post_obj