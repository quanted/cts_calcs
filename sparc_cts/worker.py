import logging
import requests
import json
import redis
import os

from sparc_calculator import SparcCalc

from cts_calcs.test_cts import views as test_views
from cts_calcs.measured_cts import views as measured_views
# from cts_app.cts_calcs.test_cts import views as test_views
# from cts_app.cts_calcs.measured_cts import views as measured_views

try:
    from cts_app.cts_calcs.smilesfilter import parseSmilesByCalculator
except ImportError as e:
    from cts_calcs.smilesfilter import parseSmilesByCalculator


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
    workflow = request.POST.get('workflow')

    # logging.info("Incoming data to SPARC: {}, {}, {} (calc, props, chemical)".format(calc, props, structure))

    post_data = {
        "calc": "sparc",
        "props": props,
        'node': node
    }

    ############### Filtered SMILES stuff!!! ###########################
    filtered_smiles = parseSmilesByCalculator(structure, "sparc") # call smilesfilter
    ###########################################################################################

    # Get melting point for sparc calculations.
    # Try Measured, then TEST..although it'll be slow
    melting_point = getMass(structure, sessionid)
    logging.warning("Using melting point: {} for SPARC calculation".format(melting_point))

    # calcObj = SparcCalc(structure, meltingpoint=melting_point)
    calcObj = SparcCalc(filtered_smiles, meltingpoint=melting_point)

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
        'request_post': request.POST,
        'run_type': run_type,
        'workflow': workflow
    }

    if run_type == 'rest':
        props = [prop]


    try:
        if 'ion_con' in props:
            response = calcObj.makeCallForPka() # response as d ict returned..
            pka_data = calcObj.getPkaResults(response)
            sparc_response.update({'data': pka_data, 'prop': 'ion_con'})
            logging.info("ion_con response: {}".format(sparc_response))
            result_json = json.dumps(sparc_response)
            redis_conn.publish(sessionid, result_json)

        if 'kow_wph' in props:
            ph = request.POST.get('ph') # get PH for logD calculation..
            response = calcObj.makeCallForLogD() # response as dict returned..
            sparc_response.update({'data': calcObj.getLogDForPH(response, ph), 'prop': 'kow_wph'})
            logging.info("kow_wph response: {}".format(sparc_response))
            result_json = json.dumps(sparc_response)
            redis_conn.publish(sessionid, result_json)

        multi_response = calcObj.makeDataRequest()
        logging.info("MULTI RESPONSE: {}".format(multi_response))
        if 'calculationResults' in multi_response:
            multi_response = calcObj.parseMultiPropResponse(multi_response['calculationResults'])
            logging.info("Parsed Multi Response: {}".format(multi_response))
            for prop_obj in multi_response:
                logging.info("PROP OBJECT: {}".format(prop_obj))
                logging.info("requested props: {}".format(props))
                logging.info("prop in props: {}".format(prop_obj['prop'] in props))
                logging.info("prop: {}".format(prop_obj['prop']))
                if prop_obj['prop'] in props and prop_obj['prop'] != 'ion_con' and prop_obj['prop'] != 'kow_wph':
                    prop_obj.update({'node': node, 'chemical': structure, 'request_post': request.POST, 'workflow': workflow, 'run_type': run_type})
                    # prop_obj.update(sparc_response)
                    logging.info("multiprop response: {}".format(prop_obj))
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


def getMass(structure, sessionid):
    """
    Gets mass of structure from Measured, tries
    TEST if not available in Measured. Returns 0.0
    if neither have mp value.
    """
    melting_point_request = {
        'calc': "measured",  # should prob be measured
        'props': ['melting_point'],
        'chemical': structure,
        'sessionid': sessionid
    }
    # todo: catch measured errors, then try epi melting point..
    request = NotDjangoRequest(melting_point_request)
    melting_point_response = measured_views.request_manager(request)

    # # convert to python dict
    try:
        melting_point = json.loads(melting_point_response.content)['data']
    except Exception as e:
        logging.warning("Error in sparc_cts/worker.py: {}".format(e))
        melting_point = 0.0

    logging.warning("MELTING POINT RESPONSE: {}".format(melting_point_response))
    logging.warning("MELTING POINT RESPONSE TYPE: {}".format(type(melting_point_response)))

    if not isinstance(melting_point, float):
        logging.warning("Trying to get MP from TEST..")
        try:
            melting_point_request['calc'] = 'test'
            request = NotDjangoRequest(melting_point_request)
            test_melting_point_response = test_views.request_manager(request)
            logging.warning("TEST MP RESPONSE CONTENT: {}".format(test_melting_point_response.content))
            melting_point = json.loads(test_melting_point_response.content)[0]['data']
            logging.warning("TEST MP VALUE: {}".format(melting_point))
        except Exception as e:
            logging.warning("Error in sparc_cts/worker.py: {}".format(e))
            melting_point = 0.0

        logging.warning("TEST MP TYPE: {}:".format(type(melting_point)))

        if not isinstance(melting_point, float):
            melting_point = 0.0
    # else:
    #     melting_point = melting_point_obj['data']

    logging.warning("MELTING POINT VALUE: {}".format(melting_point))

    return melting_point
    

class NotDjangoRequest(object):
    """
    patch for freeing celery from django while calc views
    are still relying on django.http Request...
    """
    def __init__(self, post_obj):
        self.POST = post_obj


