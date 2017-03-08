# from __future__ import absolute_import
from django.http import HttpResponse, HttpRequest
import requests
import jchem_rest  # __future__ import importerror for this...
# from cts_app.cts_calcs.chemaxon_cts import jchem_rest
import logging
import json
import redis
import jchem_calculator as JchemProperty
# from cts_app.cts_calcs.chemaxon_cts.jchem_calculator import JchemProperty
import os

logging.warning("chemaxon worker importing data_walks...")

from cts_app import cts_calcs
logging.warning("cts_calcs dir: {}".format(dir(cts_calcs)))

# from .cts_api import cts_calcs

# potential fix for cts_celery and cts_app using cts_calcs and how they import
try:
    logging.warning("trying cts_api.cts_calcs data_walks import...")
    from cts_api.cts_calcs import data_walks
except ImportError as e:
    logging.warning("trying cts_calcs data_walks import...")
    from cts_calcs import data_walks


methods = ['KLOP', 'VG', 'PHYS']

redis_hostname = os.environ.get('REDIS_HOSTNAME')
logging.warning("CHEMAXON VIEWS REDIS HOSTNAME: {}".format(redis_hostname))

redis_conn = redis.StrictRedis(host=redis_hostname, port=6379, db=0)


def request_manager(request):
    """
    Redirects request from frontend to the
    approriate jchem web service. All calls going
    from browser to cts for jchem go through
    here first.
    Expects: data for requested service, and
    which name of service to call
    Format: {"service": "", "data": {usual POST data}}
    """


    service = request.POST.get('service')
    chemical = request.POST.get('chemical')
    prop = request.POST.get('prop')
    sessionid = request.POST.get('sessionid')  # None if it doesn't exist
    method = request.POST.get('method')
    ph = request.POST.get('ph')
    node = request.POST.get('node') # TODO: keep track of nodes without bringing node obj along for the ride!
    calc = request.POST.get('calc')
    run_type = request.POST.get('run_type')
    workflow = request.POST.get('workflow')
    request_post = request.POST
    mass = request.POST.get('mass')

    logging.warning("REQUEST: {}".format(request.POST))

    try:
        props = request.POST.get("props[]")
        if not props:
            props = request.POST.getlist("props")
    except AttributeError:
        props = request.POST.get("props")


    # session = FuturesSession()  # currently not used..
    session = None


    if service == 'getTransProducts':
        # getTransProducts chemaxon service via ws..
        request = {
            'structure': chemical,
            'generationLimit': 1,  # make sure to get this from front end
            'populationLimit': 0,
            'likelyLimit': 0.001,
            'transformationLibraries': ['human_biotransformation'],  # TODO: Use UI lib choices
            'excludeCondition': ""  # 'generateImages': False
        }

        response = jchem_rest.getTransProducts(request)
        data_walks.j = 0
        data_walks.metID = 0
        results = data_walks.recursive(response, 1)

        data_obj = {
            'calc': "chemaxon", 
            'prop': "products",
            'node': node,
            'data': json.loads(results),
            'chemical': chemical,
            'workflow': 'gentrans',
            'run_type': 'batch'
        }
        
        result_json = json.dumps(data_obj)
        logging.info("publishing to redis")
        redis_conn.publish(sessionid, result_json)

    elif service == 'getSpeciationData':

        # TODO: Get this through jchem_rest and not the chemspec model!!!

        # speciation data from chemspec model, for batch via ws
        # from cts_app.models.chemspec import chemspec_output

        spec_inputs = request.POST.get('speciation_inputs')

        model_params = {
            'run_type': 'single',
            # 'run_type': 'rest',
            # 'run_type': 'batch',
            'chem_struct': None,
            'smiles': chemical,
            'orig_smiles': None,
            'iupac': None,
            'formula': None,
            'mass': None,
            'get_pka': True,
            'get_taut': None,
            'get_stereo': None,
            'pKa_decimals': None,
            'pKa_pH_lower': None,
            'pKa_pH_upper': None,
            'pKa_pH_increment': None,
            'pH_microspecies': None,
            'isoelectricPoint_pH_increment': None,
            'tautomer_maxNoOfStructures': None,
            'tautomer_pH': None,
            'stereoisomers_maxNoOfStructures': None
        }

        for key, value in model_params.items():
            if key in spec_inputs:
                model_params.update({key: spec_inputs[key]})


        # request = HttpRequest()
        # request.POST = model_params
        # request.method = "POST"  # required because of chemspec_output @request_POST???


        # TODO: CTS API URL has env var somewhere...
        # chemspec_obj = chemspec_output.chemspecOutputPage(request)
        # speciation_response = requests.post('http://localhost:8000/cts/rest/speciation/run', data=json.dumps(model_params))
        speciation_response = requests.post(
                                os.environ.get('CTS_REST_SERVER') + '/cts/rest/speciation/run', 
                                data=json.dumps(model_params), 
                                allow_redirects=True,
                                verify=False)

        speciation_data = json.loads(speciation_response.content)

        logging.warning("SPECIATION DATA RECEIVED: {}".format(speciation_data))

        data_obj = {
            'calc': "chemaxon", 
            'prop': "speciation_results",
            'node': node,
            # 'data': chemspec_obj.jchemDictResults,
            # 'data': speciation_data,
            'chemical': chemical,
            'workflow': 'chemaxon',
            'run_type': 'batch'
        }
        data_obj.update(speciation_data)

        result_json = json.dumps(data_obj)
        redis_conn.publish(sessionid, result_json)

    else:
        getPchemPropData(chemical, sessionid, method, ph, node, calc, run_type, workflow, props, session, request_post, mass)



def getPchemPropData(chemical, sessionid, method, ph, node, calc, run_type, workflow, props, session, request_post=None, mass=None):

    postData = {
        'calc': calc,
        # 'props': props
    }

    logging.warning("post data: {}".format(postData))

    if props:
        chemaxon_results = []
        for prop in props:

            logging.info("requesting chemaxon {} data".format(prop))

            data_obj = {
                'calc': calc,
                'prop':prop,
                'node': node,
                'chemical': chemical,
                'run_type': run_type,
                'workflow': workflow
            }

            # if run_type:
            #     data_obj.update({'run_type': run_type})

            try:

                data_obj.update({'request_post': request_post})

                if prop == 'kow_wph' or prop == 'kow_no_ph':
                    for method in methods:
                        
                        new_data_obj = {}
                        new_data_obj.update({
                            'calc': calc,
                            'prop': prop,
                            'chemical': chemical,
                            'node': node,
                            'request_post': request_post,
                            'run_type': run_type,
                            'workflow': workflow
                        })

                        results = getJchemPropData(chemical, prop, ph, method, sessionid, node, session)
                        new_data_obj.update({'data': results['data'], 'method': method})

                        logging.info("chemaxon results: {}".format(results))

                        result_json = json.dumps(new_data_obj)
                        redis_conn.publish(sessionid, result_json)

                else:
                    results = getJchemPropData(chemical, prop, ph, None, sessionid, node, session, mass)
                    data_obj.update({'data': results['data']})

                    logging.info("chemaxon results: {}".format(results))

                    result_json = json.dumps(data_obj)
                    redis_conn.publish(sessionid, result_json)

            except Exception as err:
                logging.warning("Exception occurred getting chemaxon data: {}".format(err))

                data_obj.update({
                    'error': "cannot reach chemaxon calculator"
                })

                redis_conn.publish(sessionid, json.dumps(data_obj))


        redis_conn.delete(sessionid)  # clear user's job cache from redis

    return


def getJchemPropData(chemical, prop, ph=7.0, method=None, sessionid=None, node=None, session=None, mass=None):
    """
    Calls jchem web services from chemaxon and
    wraps data in a CTS data object (keys: calc, prop, method, data)
    """

    resultDict = {"calc": "chemaxon", "prop": prop}

    result = ""
    if prop == 'water_sol':
        propObj = JchemProperty.getPropObject('solubility')
        propObj.makeDataRequest(chemical, None, session)
        result = propObj.getIntrinsicSolubility()
    elif prop == 'ion_con':
        propObj = JchemProperty.getPropObject('pKa')
        propObj.makeDataRequest(chemical, None, session)

        pkas = propObj.getMostAcidicPka() + propObj.getMostBasicPka()
        pkas.sort()

        result = {'pKa': pkas}
        # result = {'pKa': propObj.getMostAcidicPka(), 'pKb': propObj.getMostBasicPka()}
    elif prop == 'kow_no_ph':
        propObj = JchemProperty.getPropObject('logP')
        propObj.makeDataRequest(chemical, method, session)
        result = propObj.getLogP()
    elif prop == 'kow_wph':
        propObj = JchemProperty.getPropObject('logD')
        propObj.makeDataRequest(chemical, method, session)
        result = propObj.getLogD(ph)
    elif prop == 'water_sol_ph':
        propObj = JchemProperty.getPropObject('solubility')
        propObj.makeDataRequest(chemical, method, session)
        result = propObj.getPHDependentSolubility(ph)
        result = propObj.convertLogToMGPERL(result, mass)
    else:
        result = None

    # ADD METHOD KEY:VALUE IF LOGD OR LOGP...
    resultDict['data'] = result
    if method:
        resultDict['method'] = method

    return resultDict