import json
from .functions import request, get_json



class Assay(object):

    @classmethod
    def from_aid(cls, aid):
        """Retrieve the Assay record for the specified AID.

        :param int aid: The PubChem Assay Identifier (AID).
        """
        record = json.loads(request(aid, 'aid', 'assay', 'description').read().decode())['PC_AssayContainer'][0]
        return cls(record)

    def __init__(self, record):
        self.record = record
        """A dictionary containing the full Assay record that all other properties are obtained from."""

    def __repr__(self):
        return 'Assay(%s)' % self.aid if self.aid else 'Assay()'

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.record == other.record

    def to_dict(self, properties=None):
        """Return a dictionary containing Assay data.

        If the properties parameter is not specified, everything is included.

        :param properties: (optional) A list of the desired properties.
        """
        if not properties:
            properties = [p for p in dir(Assay) if isinstance(getattr(Assay, p), property)]
        return {p: getattr(self, p) for p in properties}

    @property
    def aid(self):
        """The PubChem Substance Idenfitier (SID)."""
        return self.record['assay']['descr']['aid']['id']

    @property
    def name(self):
        """The short assay name, used for display purposes."""
        return self.record['assay']['descr']['name']

    @property
    def description(self):
        """Description"""
        return self.record['assay']['descr']['description']

    @property
    def project_category(self):
        """A category to distinguish projects funded through MLSCN, MLPCN or from literature.

        Possible values include mlscn, mlpcn, mlscn-ap, mlpcn-ap, literature-extracted, literature-author,
        literature-publisher, rnaigi.
        """
        if 'project_category' in self.record['assay']['descr']:
            return self.record['assay']['descr']['project_category']

    @property
    def comments(self):
        """Comments and additional information."""
        return [comment for comment in self.record['assay']['descr']['comment'] if comment]

    @property
    def results(self):
        """A list of dictionaries containing details of the results from this Assay."""
        return self.record['assay']['descr']['results']

    @property
    def target(self):
        """A list of dictionaries containing details of the Assay targets."""
        if 'target' in self.record['assay']['descr']:
            return self.record['assay']['descr']['target']

    @property
    def revision(self):
        """Revision identifier for textual description."""
        return self.record['assay']['descr']['revision']

    @property
    def aid_version(self):
        """Incremented when the original depositor updates the record."""
        return self.record['assay']['descr']['aid']['version']


def get_assays(identifier, namespace='aid', **kwargs):
    """Retrieve the specified assay records from PubChem.

    :param identifier: The assay identifier to use as a search query.
    :param namespace: (optional) The identifier type.
    """
    results = get_json(identifier, namespace, 'assay', 'description', **kwargs)
    return [Assay(r) for r in results['PC_AssayContainer']] if results else []
