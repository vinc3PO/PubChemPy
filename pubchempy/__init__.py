import logging
from .functions import (get_json, get_sdf, get_sids, get_properties,
                        get_synonyms, get_aids, get_cids, get_all_sources)
from .compound import Compound, get_compounds
from .substance import Substance, get_substances
from .assay import Assay, get_assays

log = logging.getLogger('pubchempy')
log.addHandler(logging.NullHandler())

__author__ = 'Matt Swain'
__email__ = 'm.swain@me.com'
__version__ = '1.0.4'
__license__ = 'MIT'