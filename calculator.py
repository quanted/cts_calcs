__author__ = 'np'

from django.template import Context, Template
import requests
import json
import logging
import os
import redis
import datetime
import pytz

# from cts_app.cts_calcs import jchem_properties
# import jchem_properties


class Calculator(object):
    """
    Skeleton class for calculators
    """

    def __init__(self, calc=None):
        """
        calc -- p-chem calculator name, 'getTransProducts', or 'getSpeciationData'
        """
        self.name = ''
        self.propMap = {}
        self.baseUrl = None
        self.urlStruct = ''
        self.results = ''
        self.headers = {'Content-Type': 'application/json'}
        self.request_timeout = 30  # default, set unique ones in calc sub classes
        self.max_retries = 3
        # self.request_type = ''  # http or ws (https or wss [todo])

        self.image_scale = 50

        self.redis_hostname = os.environ.get('REDIS_HOSTNAME')
        self.redis_port = os.environ.get('REDIS_PORT')
        self.redis_conn = redis.StrictRedis(host=self.redis_hostname, port=self.redis_port, db=0)

        self.jchem_server_url = os.environ.get('CTS_JCHEM_SERVER', 'localhost:8080')
        self.efs_server_url = os.environ.get('CTS_EFS_SERVER', 'localhost:8080')

        # jchem ws urls:
        self.export_endpoint = '/webservices/rest-v0/util/calculate/molExport'
        self.detail_endpoint = '/webservices/rest-v0/util/detail'
        # self.hydro_endpoint = '/webservices/rest-v0/util/convert/hydrogenizer'
        # self.standardizer_endpoint = '/webservices/rest-v0/util/convert/standardizer'

        # CTSWS (formerlly EFS) metabolizer endpoints:
        self.efs_metabolizer_endpoint = '/ctsws/rest/metabolizer'
        self.efs_standardizer_endpoint = '/ctsws/rest/standardizer'

        # cts p-chem properties
        self.pchem_props = [
            'boiling_point',
            'melting_point',
            'water_sol',
            'vapor_press',
            'mol_diss',
            'ion_con',
            'henrys_law_con',
            'kow_no_ph',
            'kow_wph',
            'kow_ph',
            'kow'
        ]

        # cts chemical information dict
        self.chemical_information = {
            'chemical': None,  # user-entered chemical (as-entered or drawn)
            'orig_smiles': None,  # original conversion to SMILES
            'smiles': None,  # SMILES after filtering, used for calculations
            'formula': None,
            'iupac': None,
            'mass': None,
            'structureData': None,  # drawn chemical structure format for MarvinSketch
            'exactMass': None,
        }

        # cts chemical information request
        self.chemical_information_request = {
            'chemical': None,
            'get_structure_data': False,
        }

        # cts api data object for p-chem data request
        self.data_obj = {
            'calc': None,
            'prop': None,
            'data': None,
            'chemical': None,
        }

        # cts p-chem request object with default key:vals.
        # can handle list of props (ws) or single prop (cts api)
        self.pchem_request = {
            'service': None,
            'chemical': None,
            'prop': None,
            'sessionid': None,
            'method': None,
            'ph': 7.0,
            'node': None,
            'calc': None,
            'run_type': None,
            'workflow': None,
            'mass': None,
            'props': [],
        }

        # cts p-chem response object with defaults, returns one prop per reponse
        self.pchem_response = {
            'chemical': None,
            'calc': None,
            'prop': None,
            'method': None,
            'run_type': None,
            'workflow': None,
            'node': None,
            'request_post': None,
            'data': None,
            'error': False,  # ehh?
        }


    def getUrl(self, prop):
        if prop in self.propMap:
            calcProp = self.propMap[prop]['urlKey']
            return self.urlStruct.format(calcProp)
        else:
            return "Error: url key not found"

    def getPropKey(self, prop):
        if prop in self.propMap:
            return self.propMap[prop]['propKey']
        else:
            return "Error: prop not found"

    def getResultKey(self, prop):
        if prop in self.propMap:
            return self.propMap[prop]['resultKey']
        else:
            return "Error: result key not found"






    ################### FUNCTION FROM CTS_REST USED BY BASICALLY EVERYTHING ###########

    def gen_jid(self):
        ts = datetime.datetime.now(pytz.UTC)
        localDatetime = ts.astimezone(pytz.timezone('US/Eastern'))
        jid = localDatetime.strftime('%Y%m%d%H%M%S%f')
        return jid





    ################ JCHEM REST STUFF (WHERE SHOULD IT GO??? CTS_REST???) ###################

    def getChemDetails(self, request_obj):
        """
        getChemDetails

        Inputs:
        chem - chemical name (format: iupac, smiles, or formula)
        Returns:
        The iupac, formula, mass, and smiles string of the chemical
        along with the mrv of the chemical (to display in marvinjs)
        """
        chemical = request_obj.get('chemical')

        chemDeatsDict = {
            "structures": [
                {"structure": chemical}
            ],
            "display": {
                "include": [
                    "structureData"
                ],
                "additionalFields": {
                    "formula": "chemicalTerms(formula)",
                    "iupac": "chemicalTerms(name)",
                    "mass": "chemicalTerms(mass)",
                    "exactMass": "chemicalTerms(exactMass)",
                    "smiles": "chemicalTerms(molString('smiles'))",
                    "cas": "chemicalTerms(molString('name:cas#'))",
                },
                "parameters": {
                    "structureData": "mrv"
                }
            }
        }

        url = self.jchem_server_url + self.detail_endpoint
        return self.web_call(url, chemDeatsDict)


    def smilesToImage(self, request_obj):
        """
        smilesToImage

        Returns image (.png) url for a 
        given SMILES
        """
        smiles = request_obj.get('smiles')
        imgScale = request_obj.get('scale')
        imgWidth = request_obj.get('width')
        imgHeight = request_obj.get('height')
        imgType = request_obj.get('type')

        # NOTE: Requesting image without width or height but scale
        # returns an image that just fits the molecule with nice resoltuion.
        # Providing width and height without scale might return a higher
        # resolution image for metabolites!
        
        request = {
            "structures": [
                {"structure": smiles}
            ],
            "display": {
                "include": ["image"],
                "parameters": {
                    "image": {
                    }
                }
            }
        }

        if imgType:
            request['display']['parameters']['image'].update({'type': imgType})
        else:
            request['display']['parameters']['image'].update({'type': 'png'})

        if imgHeight:
            # these are metabolites in the space tree:
            request['display']['parameters']['image'].update({"width": imgWidth, "height": imgHeight})
        else:
            request['display']['parameters']['image'].update({'width': imgWidth, 'scale': imgScale})

        url = self.jchem_server_url + self.detail_endpoint
        imgData = self.web_call(url, request)  # get response from jchem ws
        return imgData  # return dict of image data


    def convertToSMILES(self, request_obj):
        """
        convertToSMILES

        Inputs: chemical as mrv, smiles, etc. (chemaxon recognized)
        Returns: SMILES string of chemical
        """
        chemStruct = request_obj.get('chemical')  # chemical in <cml> format (marvin sketch)
        data = {
            "structure": chemStruct,
            # "inputFormat": "mrv",
            "parameters": "smiles"
        }
        url = self.jchem_server_url + self.export_endpoint
        return self.web_call(url, data)  # get responset))


    # def getTransProducts(self, request_obj):
    #     """
    #     Makes request to metabolizer
    #     """
    #     url = self.efs_server_url + self.efs_metabolizer_endpoint
    #     return self.web_call(url, request_obj)


    def getStructInfo(self, structure):
        """
        Appends structure info to image url
        Input: structure in .mrv format
        Output: dict with structure's info (i.e., formula, iupac, mass, smiles),
        or dict with aforementioned keys but None values
        """
        structDict = self.getChemDetails({"chemical": structure, "addH": True})
        infoDictKeys = ['formula', 'iupac', 'mass', 'smiles','exactMass']
        infoDict = {key: None for key in infoDictKeys}  # init dict with infoDictKeys and None vals
        struct_root = {}  # root of data in structInfo
        if 'data' in structDict:
            struct_root = structDict['data'][0]
            infoDict.update({
                "formula": struct_root['formula'],
                "iupac": struct_root['iupac'],
                "mass": struct_root['mass'],
                "smiles": struct_root['smiles'],
                'exactMass': struct_root['exactMass']
            })
        return infoDict


    def getMass(self, request_obj):
        """
        get mass of structure from jchem ws
        """
        chemical = request_obj.get('chemical')
        logging.info("jchem_rest getting mass for {}".format(chemical))
        post_data = {
            "structures": [
                {"structure": chemical}
            ],
            "display": {
                "include": [
                    "structureData"
                ],
                "additionalFields": {
                    "mass": "chemicalTerms(mass)"
                },
                "parameters": {
                    "structureData": "smiles"
                }
            }
        }
        url = self.jchem_server_url + self.detail_endpoint
        return self.web_call(url, post_data)


    def web_call(self, url, data):
        """
        Makes the request to a specified URL
        and POST data. Returns resonse data as dict
        """

        # TODO: Deal with errors more granularly... 403, 500, etc.
        try:
            if data == None:
                response = requests.get(url, timeout=self.request_timeout)
            else:
                response = requests.post(url, data=json.dumps(data), headers=self.headers, timeout=self.request_timeout)
            return json.loads(response.content)
        except requests.exceptions.RequestException as e:
            logging.warning("error at web call: {} /error".format(e))
            raise e





    ################# TRANSFORMATION PRODUCTS STUFF (DATA_WALKS) #######################

    def nodeWrapper(self, smiles, height, width, scale, key=None, img_type=None, isProduct=None):
        """
        Wraps image html tag around
        the molecule's image source
        Inputs: smiles, height, width, scale, key
        Returns: html of wrapped image
        """

        # 1. Get image from smiles
        post = {
            "smiles": smiles,
            "scale": scale,
            "height": height,
            "width": width,
            # "type": img_type
        }

        if img_type:
            post.update({'type': img_type})

        results = self.smilesToImage(post)

        # 2. Get imageUrl out of results
        img, imgScale = '', ''
        if 'data' in results:
            root = results['data'][0]['image']
            if 'image' in root:
                img = root['image']

        # 3. Wrap imageUrl with <img>
        # <img> wrapper for image byte string:
        if img_type and img_type == 'svg':

            if key:
                html = '<div style="background-color:white;"' + 'id="' + str(key) + '">' + img + '</div>'            
            else:
                html = '<div style="background-color:white;">' + img + '</div>'
            
        else:
            html = self.imgTmpl(isProduct).render(Context(dict(smiles=smiles, img=img, height=height, width=width, scale=scale, key=key)))

        # wrapper for SVG type images:
        # html = data['data'][0]['image']['image']  # expecting svg type now (could affect popups)

        # do the <svg> elements need class and id like imgTmpl below???

        return html


    def imgTmpl(self, isProduct):
        if isProduct:
            imgTmpl = """
            <img class="metabolite" id="{{key|default:""}}"
                alt="{{smiles}}" src="data:image/png;base64,{{img}}"
                width={{width}} height={{height}} /> 
            """        
        else:
            imgTmpl = """
            <img class="metabolite" id="{{key|default:""}}"
                alt="{{smiles}}" src="data:image/png;base64,{{img}}"
                width={{width}} height={{height}} hidden /> 
            """
        return Template(imgTmpl)


    def popupBuilder(self, root, paramKeys, molKey=None, header=None, isProduct=False):
        """
        Wraps molecule data (e.g., formula, iupac, mass, 
        smiles, image) for hover-over popups in chemspec
        and gentrans outputs.

        Inputs:
        root - dictionary of items to wrap in table
        paramKeys - keys to use for building table
        molKey - (optional) add id to wrap table
        header - (optional) add header above key/values 

        Returns: dictionary where html key is 
        the wrapped html and the other keys are
        same as the input keys
        """

        # propKeys = ['smiles', 'accumulation', 'production', 'transmissivity', 'generation']
        dataProps = {key: None for key in paramKeys}  # metabolite properties

        html = '<div id="{}_div" class="nodeWrapDiv"><div class="metabolite_img" style="float:left;">'.format(molKey)
        # html += nodeWrapper(root['smiles'], None, 250, 150)

        # smiles, height, width, scale, key=None, img_type=None

        if isProduct:
            html += self.nodeWrapper(root['smiles'], None, 250, self.image_scale, molKey, 'png')  # hidden png for pdf    
        else:
            # html += nodeWrapper(root['smiles'], None, 250, image_scale)  # Molecular Info image, metabolites output
            html += self.nodeWrapper(root['smiles'], None, 250, self.image_scale, molKey, 'svg')  # svg popups for chemspec and gentrans outputs
            html += self.nodeWrapper(root['smiles'], None, 250, self.image_scale, molKey, None)  # hidden png for pdf

        html += '</div>'

        if molKey:
            html += '<table class="ctsTableStylin" id="{}_table">'.format(molKey)
        else:
            html += '<table class="ctsTableStylin">'

        if header:
            html += '<tr class="header"><th colspan="2">' + header + '</th></tr>'

        for key, value in root.items():
            if key in paramKeys:

                # Convert other types (e.g., float, int) to string
                if not isinstance(value, unicode) and not (isinstance(value, str)):
                    
                    if key == 'exactMass':
                        value = str(value)
                    else:
                        value = str(round(float(value), 3))
                # value = str(value)

                dataProps[key] = value

                html += '<tr><td>' + key + '</td>'
                html += '<td>' + value + '</td></tr>'
        html += '</table></div>'

        dataProps["html"] = html

        return dataProps