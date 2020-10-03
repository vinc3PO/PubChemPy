
# -*- coding: utf-8 -*-
"""
PubChemPy

Python interface for the PubChem PUG REST service.
https://github.com/mcs07/PubChemPy
"""

import json
import os
import time

from .mapper import PROPERTY_MAP

from urllib.error import HTTPError
from urllib.parse import quote, urlencode
from urllib.request import urlopen
from .errors import PubChemHTTPError, NotFoundError
import re
import logging



API_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
API_VIEW = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound'



text_types = str, bytes


def request(identifier, namespace='cid', domain='compound', operation=None, output='JSON', searchtype=None, **kwargs):
    """
    Construct API request from parameters and return the response.

    Full specification at http://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    """
    if not identifier:
        raise ValueError('identifier/cid cannot be None')
    # If identifier is a list, join with commas into string
    if isinstance(identifier, int):
        identifier = str(identifier)
    if not isinstance(identifier, text_types):
        identifier = ','.join(str(x) for x in identifier)
    # Filter None values from kwargs
    kwargs = dict((k, v) for k, v in kwargs.items() if v is not None)
    # Build API URL
    urlid, postdata = None, None
    if namespace == 'sourceid':
        identifier = identifier.replace('/', '.')
    if namespace in ['listkey', 'formula', 'sourceid'] \
            or searchtype == 'xref' \
            or (searchtype and namespace == 'cid') or domain == 'sources':
        urlid = quote(identifier.encode('utf8'))
    else:
        postdata = urlencode([(namespace, identifier)]).encode('utf8')
    comps = filter(None, [API_BASE, domain, searchtype, namespace, urlid, operation, output])
    apiurl = '/'.join(comps)
    if kwargs:
        apiurl += '?%s' % urlencode(kwargs)
    # Make request
    try:
        # log.debug('Request URL: %s', apiurl)
        # log.debug('Request data: %s', postdata)
        response = urlopen(apiurl, postdata)
        return response
    except HTTPError as e:
        raise PubChemHTTPError(e)


def get(identifier, namespace='cid', domain='compound', operation=None, output='JSON', searchtype=None, **kwargs):
    """Request wrapper that automatically handles async requests."""
    if (searchtype and searchtype != 'xref') or namespace in ['formula']:
        response = request(identifier, namespace, domain, None, 'JSON', searchtype, **kwargs).read()
        status = json.loads(response.decode())
        if 'Waiting' in status and 'ListKey' in status['Waiting']:
            identifier = status['Waiting']['ListKey']
            namespace = 'listkey'
            while 'Waiting' in status and 'ListKey' in status['Waiting']:
                time.sleep(2)
                response = request(identifier, namespace, domain, operation, 'JSON', **kwargs).read()
                status = json.loads(response.decode())
            if not output == 'JSON':
                response = request(identifier, namespace, domain, operation, output, searchtype, **kwargs).read()
    else:
        response = request(identifier, namespace, domain, operation, output, searchtype, **kwargs).read()
    return response


def get_json(identifier, namespace='cid', domain='compound', operation=None, searchtype=None, **kwargs):
    """Request wrapper that automatically parses JSON response and supresses NotFoundError."""
    try:
        return json.loads(get(identifier, namespace, domain, operation, 'JSON', searchtype, **kwargs).decode())
    except NotFoundError as e:
        # log.info(e)
        return None

def get_sdf(identifier, namespace='cid', domain='compound',operation=None, searchtype=None, **kwargs):
    """Request wrapper that automatically parses SDF response and supresses NotFoundError."""
    try:
        return get(identifier, namespace, domain, operation, 'SDF', searchtype, **kwargs).decode()
    except NotFoundError as e:
        # log.info(e)
        return None


def get_properties(properties, identifier, namespace='cid', searchtype=None, as_dataframe=False, **kwargs):
    """Retrieve the specified properties from PubChem.

    :param identifier: The compound, substance or assay identifier to use as a search query.
    :param namespace: (optional) The identifier type.
    :param searchtype: (optional) The advanced search type, one of substructure, superstructure or similarity.
    :param as_dataframe: (optional) Automatically extract the properties into a pandas :class:`~pandas.DataFrame`.
    """
    if isinstance(properties, text_types):
        properties = properties.split(',')
    properties = ','.join([PROPERTY_MAP.get(p, p) for p in properties])
    properties = 'property/%s' % properties
    results = get_json(identifier, namespace, 'compound', properties, searchtype=searchtype, **kwargs)
    results = results['PropertyTable']['Properties'] if results else []
    if as_dataframe:
        import pandas as pd
        return pd.DataFrame.from_records(results, index='CID')
    return results


def get_synonyms(identifier, namespace='cid', domain='compound', searchtype=None, **kwargs):
    results = get_json(identifier, namespace, domain, 'synonyms', searchtype=searchtype, **kwargs)
    return results['InformationList']['Information'] if results else []


def get_cids(identifier, namespace='name', domain='compound', searchtype=None, **kwargs):
    results = get_json(identifier, namespace, domain, 'cids', searchtype=searchtype, **kwargs)
    if not results:
        return []
    elif 'IdentifierList' in results:
        return results['IdentifierList']['CID']
    elif 'InformationList' in results:
        return results['InformationList']['Information']


def get_sids(identifier, namespace='cid', domain='compound', searchtype=None, **kwargs):
    results = get_json(identifier, namespace, domain, 'sids', searchtype=searchtype, **kwargs)
    if not results:
        return []
    elif 'IdentifierList' in results:
        return results['IdentifierList']['SID']
    elif 'InformationList' in results:
        return results['InformationList']['Information']


def get_aids(identifier, namespace='cid', domain='compound', searchtype=None, **kwargs):
    results = get_json(identifier, namespace, domain, 'aids', searchtype=searchtype, **kwargs)
    if not results:
        return []
    elif 'IdentifierList' in results:
        return results['IdentifierList']['AID']
    elif 'InformationList' in results:
        return results['InformationList']['Information']


def get_all_sources(domain='substance'):
    """Return a list of all current depositors of substances or assays."""
    results = json.loads(get(domain, None, 'sources').decode())
    return results['InformationList']['SourceName']


def download(outformat, path, identifier, namespace='cid', domain='compound', operation=None, searchtype=None,
             overwrite=False, **kwargs):
    """Format can be  XML, ASNT/B, JSON, SDF, CSV, PNG, TXT."""
    response = get(identifier, namespace, domain, operation, outformat, searchtype, **kwargs)
    if not overwrite and os.path.isfile(path):
        raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % path)
    with open(path, 'wb') as f:
        f.write(response)


def _parse_prop(search, proplist):
    """Extract property value from record using the given urn search filter."""
    props = [i for i in proplist if all(item in i['urn'].items() for item in search.items())]
    if len(props) > 0:
        return props[0]['value'][list(props[0]['value'].keys())[0]]


def request_SDS(cid):
    if not cid:
        raise ValueError('identifier/cid cannot be None')
    # Make request
    try:
        # log.debug('Request URL: %s', apiurl)
        # log.debug('Request data: %s', postdata)
        response = urlopen(API_VIEW + '/{}/JSON?heading=safety+and+hazards'.format(cid))
        return _parse_sds(json.loads(response.read().decode()))
    except HTTPError as e:
        raise PubChemHTTPError(e)


def _parse_sds(result):
    info = result["Record"]["Section"][0]["Section"][0]['Section'][0]['Information']
    pictogram = []
    hazard = []
    precautionary = []
    for resp in info:
        if resp['Name'] == 'Pictogram(s)':
            for picto in resp['Value']['StringWithMarkup'][0]['Markup']:
                pictObject = {"icon": picto['URL'].split('/')[-1], "string": picto['Extra']}
                if pictObject not in pictogram:
                    pictogram.append(pictObject)
        elif resp['Name'] == 'GHS Hazard Statements':
            for stat in resp["Value"]['StringWithMarkup']:
                if stat['String'][0] == "H" and stat['String'][1:4].isnumeric():
                    if stat['String'][0:4] not in hazard:
                        hazard.append(stat['String'][0:4])
                        hazard.sort()
        elif resp['Name'] == 'Precautionary Statement Codes':
            for codes in resp["Value"]["StringWithMarkup"]:
                r1 = re.findall(r"[P0-9+]{3,30}", codes['String'])
                if r1:
                    for ps in r1:
                        if ps not in precautionary:
                            precautionary.append(ps)
                            precautionary.sort()
    return pictogram, hazard, precautionary

