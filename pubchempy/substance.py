import json
from .functions import get_json, request
from .mapper import CompoundIdType
from .decorators import memoized_property
from .compound import Compound
from .logger import createLogger

log = createLogger(__name__)

class Substance(object):
    """Corresponds to a single record from the PubChem Substance database.

    The PubChem Substance database contains all chemical records deposited in PubChem in their most raw form, before
    any significant processing is applied. As a result, it contains duplicates, mixtures, and some records that don't
    make chemical sense. This means that Substance records contain fewer calculated properties, however they do have
    additional information about the original source that deposited the record.

    The PubChem Compound database is constructed from the Substance database using a standardization and deduplication
    process. Hence each Compound may be derived from a number of different Substances.
    """

    @classmethod
    def from_sid(cls, sid):
        """Retrieve the Substance record for the specified SID.

        :param int sid: The PubChem Substance Identifier (SID).
        """
        record = json.loads(request(sid, 'sid', 'substance').read().decode())['PC_Substances'][0]
        return cls(record)

    def __init__(self, record):
        self.record = record
        """A dictionary containing the full Substance record that all other properties are obtained from."""

    def __repr__(self):
        return 'Substance(%s)' % self.sid if self.sid else 'Substance()'

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.record == other.record

    def to_dict(self, properties=None):
        """Return a dictionary containing Substance data.

        If the properties parameter is not specified, everything except cids and aids is included. This is because the
        aids and cids properties each require an extra request to retrieve.

        :param properties: (optional) A list of the desired properties.
        """
        if not properties:
            skip = {'deposited_compound', 'standardized_compound', 'cids', 'aids'}
            properties = [p for p in dir(Substance) if isinstance(getattr(Substance, p), property) and p not in skip]
        return {p: getattr(self, p) for p in properties}

    def to_series(self, properties=None):
        """Return a pandas :class:`~pandas.Series` containing Substance data.

        If the properties parameter is not specified, everything except cids and aids is included. This is because the
        aids and cids properties each require an extra request to retrieve.

        :param properties: (optional) A list of the desired properties.
        """
        import pandas as pd
        return pd.Series(self.to_dict(properties))

    @property
    def sid(self):
        """The PubChem Substance Idenfitier (SID)."""
        return self.record['sid']['id']

    @property
    def synonyms(self):
        """A ranked list of all the names associated with this Substance."""
        if 'synonyms' in self.record:
            return self.record['synonyms']

    @property
    def source_name(self):
        """The name of the PubChem depositor that was the source of this Substance."""
        return self.record['source']['db']['name']

    @property
    def source_id(self):
        """Unique ID for this Substance within those from the same PubChem depositor source."""
        return self.record['source']['db']['source_id']['str']

    @property
    def standardized_cid(self):
        """The CID of the Compound that was produced when this Substance was standardized.

        May not exist if this Substance was not standardizable.
        """
        for c in self.record['compound']:
            if c['id']['type'] == CompoundIdType.STANDARDIZED:
                return c['id']['id']['cid']

    @memoized_property
    def standardized_compound(self):
        """Return the :class:`~pubchempy.Compound` that was produced when this Substance was standardized.

        Requires an extra request. Result is cached.
        """
        for c in self.record['compound']:
            if c['id']['type'] == CompoundIdType.STANDARDIZED:
                return Compound.from_cid(c['id']['id']['cid'])

    @property
    def deposited_compound(self):
        """Return a :class:`~pubchempy.Compound` produced from the unstandardized Substance record as deposited.

        The resulting :class:`~pubchempy.Compound` will not have a ``cid`` and will be missing most properties.
        """
        for c in self.record['compound']:
            if c['id']['type'] == CompoundIdType.DEPOSITED:
                return Compound(c)

    @memoized_property
    def cids(self):
        """A list of all CIDs for Compounds that were produced when this Substance was standardized.

        Requires an extra request. Result is cached."""
        results = get_json(self.sid, 'sid', 'substance', 'cids')
        return results['InformationList']['Information'][0]['CID'] if results else []

    @memoized_property
    def aids(self):
        """A list of all AIDs for Assays associated with this Substance.

        Requires an extra request. Result is cached."""
        results = get_json(self.sid, 'sid', 'substance', 'aids')
        return results['InformationList']['Information'][0]['AID'] if results else []

def get_substances(identifier, namespace='sid', as_dataframe=False, **kwargs):
    """Retrieve the specified substance records from PubChem.

    :param identifier: The substance identifier to use as a search query.
    :param namespace: (optional) The identifier type, one of sid, name or sourceid/<source name>.
    :param as_dataframe: (optional) Automatically extract the :class:`~pubchempy.Substance` properties into a pandas
                         :class:`~pandas.DataFrame` and return that.
    """
    results = get_json(identifier, namespace, 'substance', **kwargs)
    substances = [Substance(r) for r in results['PC_Substances']] if results else []
    if as_dataframe:
        return substances_to_frame(substances)
    return substances



def substances_to_frame(substances, properties=None):
    """Construct a pandas :class:`~pandas.DataFrame` from a list of :class:`~pubchempy.Substance` objects.

    Optionally specify a list of the desired :class:`~pubchempy.Substance` properties.
    """
    import pandas as pd
    if isinstance(substances, Substance):
        substances = [substances]
    properties = set(properties) | set(['sid']) if properties else None
    return pd.DataFrame.from_records([s.to_dict(properties) for s in substances], index='sid')