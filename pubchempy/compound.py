import json
from .functions import get_json, request, _parse_prop
from .decorators import deprecated, memoized_property
from .mapper import ELEMENTS, CoordinateType, BondType
from .errors import ResponseParseError
from itertools import zip_longest
#from . import log

class Atom(object):
    """Class to represent an atom in a :class:`~pubchempy.Compound`."""

    def __init__(self, aid, number, x=None, y=None, z=None, charge=0):
        """Initialize with an atom ID, atomic number, coordinates and optional change.

        :param int aid: Atom ID
        :param int number: Atomic number
        :param float x: X coordinate.
        :param float y: Y coordinate.
        :param float z: (optional) Z coordinate.
        :param int charge: (optional) Formal charge on atom.
        """
        self.aid = aid
        """The atom ID within the owning Compound."""
        self.number = number
        """The atomic number for this atom."""
        self.x = x
        """The x coordinate for this atom."""
        self.y = y
        """The y coordinate for this atom."""
        self.z = z
        """The z coordinate for this atom. Will be ``None`` in 2D Compound records."""
        self.charge = charge
        """The formal charge on this atom."""

    def __repr__(self):
        return 'Atom(%s, %s)' % (self.aid, self.element)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.aid == other.aid and self.element == other.element and
                self.x == other.x and self.y == other.y and self.z == other.z and self.charge == other.charge)

    @deprecated('Dictionary style access to Atom attributes is deprecated')
    def __getitem__(self, prop):
        """Allow dict-style access to attributes to ease transition from when atoms were dicts."""
        if prop in {'element', 'x', 'y', 'z', 'charge'}:
            return getattr(self, prop)
        raise KeyError(prop)

    @deprecated('Dictionary style access to Atom attributes is deprecated')
    def __setitem__(self, prop, val):
        """Allow dict-style setting of attributes to ease transition from when atoms were dicts."""
        setattr(self, prop, val)

    @deprecated('Dictionary style access to Atom attributes is deprecated')
    def __contains__(self, prop):
        """Allow dict-style checking of attributes to ease transition from when atoms were dicts."""
        if prop in {'element', 'x', 'y', 'z', 'charge'}:
            return getattr(self, prop) is not None
        return False

    @property
    def element(self):
        """The element symbol for this atom."""
        return ELEMENTS.get(self.number, None)

    def to_dict(self):
        """Return a dictionary containing Atom data."""
        data = {'aid': self.aid, 'number': self.number, 'element': self.element}
        for coord in {'x', 'y', 'z'}:
            if getattr(self, coord) is not None:
                data[coord] = getattr(self, coord)
        if self.charge is not 0:
            data['charge'] = self.charge
        return data

    def set_coordinates(self, x, y, z=None):
        """Set all coordinate dimensions at once."""
        self.x = x
        self.y = y
        self.z = z

    @property
    def coordinate_type(self):
        """Whether this atom has 2D or 3D coordinates."""
        return '2d' if self.z is None else '3d'


class Bond(object):
    """Class to represent a bond between two atoms in a :class:`~pubchempy.Compound`."""

    def __init__(self, aid1, aid2, order=BondType.SINGLE, style=None):
        """Initialize with begin and end atom IDs, bond order and bond style.

        :param int aid1: Begin atom ID.
        :param int aid2: End atom ID.
        :param int order: Bond order.
        """
        self.aid1 = aid1
        """ID of the begin atom of this bond."""
        self.aid2 = aid2
        """ID of the end atom of this bond."""
        self.order = order
        """Bond order."""
        self.style = style
        """Bond style annotation."""

    def __repr__(self):
        return 'Bond(%s, %s, %s)' % (self.aid1, self.aid2, self.order)

    def __eq__(self, other):
        return (isinstance(other, type(self)) and self.aid1 == other.aid1 and self.aid2 == other.aid2 and
                self.order == other.order and self.style == other.style)

    @deprecated('Dictionary style access to Bond attributes is deprecated')
    def __getitem__(self, prop):
        """Allow dict-style access to attributes to ease transition from when bonds were dicts."""
        if prop in {'order', 'style'}:
            return getattr(self, prop)
        raise KeyError(prop)

    @deprecated('Dictionary style access to Bond attributes is deprecated')
    def __setitem__(self, prop, val):
        """Allow dict-style setting of attributes to ease transition from when bonds were dicts."""
        setattr(self, prop, val)

    @deprecated('Dictionary style access to Atom attributes is deprecated')
    def __contains__(self, prop):
        """Allow dict-style checking of attributes to ease transition from when bonds were dicts."""
        if prop in {'order', 'style'}:
            return getattr(self, prop) is not None
        return False

    @deprecated('Dictionary style access to Atom attributes is deprecated')
    def __delitem__(self, prop):
        """Delete the property prop from the wrapped object."""
        if not hasattr(self.__wrapped, prop):
            raise KeyError(prop)
        delattr(self.__wrapped, prop)

    def to_dict(self):
        """Return a dictionary containing Bond data."""
        data = {'aid1': self.aid1, 'aid2': self.aid2, 'order': self.order}
        if self.style is not None:
            data['style'] = self.style
        return data


class Compound(object):
    """Corresponds to a single record from the PubChem Compound database.

    The PubChem Compound database is constructed from the Substance database using a standardization and deduplication
    process. Each Compound is uniquely identified by a CID.
    """
    def __init__(self, record):
        """Initialize with a record dict from the PubChem PUG REST service.

        For most users, the ``from_cid()`` class method is probably a better way of creating Compounds.

        :param dict record: A compound record returned by the PubChem PUG REST service.
        """
        self._record = None
        self._atoms = {}
        self._bonds = {}
        self.record = record

    @property
    def record(self):
        """The raw compound record returned by the PubChem PUG REST service."""
        return self._record

    @record.setter
    def record(self, record):
        self._record = record
        #log.debug('Created %s' % self)
        self._setup_atoms()
        self._setup_bonds()

    def _setup_atoms(self):
        """Derive Atom objects from the record."""
        # Delete existing atoms
        self._atoms = {}
        # Create atoms
        aids = self.record['atoms']['aid']
        elements = self.record['atoms']['element']
        if not len(aids) == len(elements):
            raise ResponseParseError('Error parsing atom elements')
        for aid, element in zip(aids, elements):
            self._atoms[aid] = Atom(aid=aid, number=element)
        # Add coordinates
        if 'coords' in self.record:
            coord_ids = self.record['coords'][0]['aid']
            xs = self.record['coords'][0]['conformers'][0]['x']
            ys = self.record['coords'][0]['conformers'][0]['y']
            zs = self.record['coords'][0]['conformers'][0].get('z', [])
            if not len(coord_ids) == len(xs) == len(ys) == len(self._atoms) or (zs and not len(zs) == len(coord_ids)):
                raise ResponseParseError('Error parsing atom coordinates')
            for aid, x, y, z in zip_longest(coord_ids, xs, ys, zs):
                self._atoms[aid].set_coordinates(x, y, z)
        # Add charges
        if 'charge' in self.record['atoms']:
            for charge in self.record['atoms']['charge']:
                self._atoms[charge['aid']].charge = charge['value']

    def _setup_bonds(self):
        """Derive Bond objects from the record."""
        self._bonds = {}
        if 'bonds' not in self.record:
            return
        # Create bonds
        aid1s = self.record['bonds']['aid1']
        aid2s = self.record['bonds']['aid2']
        orders = self.record['bonds']['order']
        if not len(aid1s) == len(aid2s) == len(orders):
            raise ResponseParseError('Error parsing bonds')
        for aid1, aid2, order in zip(aid1s, aid2s, orders):
            self._bonds[frozenset((aid1, aid2))] = Bond(aid1=aid1, aid2=aid2, order=order)
        # Add styles
        if 'coords' in self.record and 'style' in self.record['coords'][0]['conformers'][0]:
            aid1s = self.record['coords'][0]['conformers'][0]['style']['aid1']
            aid2s = self.record['coords'][0]['conformers'][0]['style']['aid2']
            styles = self.record['coords'][0]['conformers'][0]['style']['annotation']
            for aid1, aid2, style in zip(aid1s, aid2s, styles):
                self._bonds[frozenset((aid1, aid2))].style = style

    @classmethod
    def from_cid(cls, cid, **kwargs):
        """Retrieve the Compound record for the specified CID.

        Usage::

            c = Compound.from_cid(6819)

        :param int cid: The PubChem Compound Identifier (CID).
        """
        record = json.loads(request(cid, **kwargs).read().decode())['PC_Compounds'][0]
        return cls(record)

    def __repr__(self):
        return 'Compound(%s)' % self.cid if self.cid else 'Compound()'

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.record == other.record

    def to_dict(self, properties=None):
        """Return a dictionary containing Compound data. Optionally specify a list of the desired properties.

        synonyms, aids and sids are not included unless explicitly specified using the properties parameter. This is
        because they each require an extra request.
        """
        if not properties:
            skip = {'aids', 'sids', 'synonyms'}
            properties = [p for p in dir(Compound) if isinstance(getattr(Compound, p), property) and p not in skip]
        return {p: [i.to_dict() for i in getattr(self, p)] if p in {'atoms', 'bonds'} else getattr(self, p) for p in properties}

    def to_series(self, properties=None):
        """Return a pandas :class:`~pandas.Series` containing Compound data. Optionally specify a list of the desired
        properties.

        synonyms, aids and sids are not included unless explicitly specified using the properties parameter. This is
        because they each require an extra request.
        """
        import pandas as pd
        return pd.Series(self.to_dict(properties))

    @property
    def cid(self):
        """The PubChem Compound Identifier (CID).

        .. note::

            When searching using a SMILES or InChI query that is not present in the PubChem Compound database, an
            automatically generated record may be returned that contains properties that have been calculated on the
            fly. These records will not have a CID property.
        """
        if 'id' in self.record and 'id' in self.record['id'] and 'cid' in self.record['id']['id']:
            return self.record['id']['id']['cid']

    @property
    def elements(self):
        """List of element symbols for atoms in this Compound."""
        return [a.element for a in self.atoms]

    @property
    def atoms(self):
        """List of :class:`Atoms <pubchempy.Atom>` in this Compound."""
        return sorted(self._atoms.values(), key=lambda x: x.aid)

    @property
    def bonds(self):
        """List of :class:`Bonds <pubchempy.Bond>` between :class:`Atoms <pubchempy.Atom>` in this Compound."""
        return sorted(self._bonds.values(), key=lambda x: (x.aid1, x.aid2))

    @memoized_property
    def synonyms(self):
        """A ranked list of all the names associated with this Compound.

        Requires an extra request. Result is cached.
        """
        if self.cid:
            results = get_json(self.cid, operation='synonyms')
            return results['InformationList']['Information'][0]['Synonym'] if results else []

    @memoized_property
    def sids(self):
        """Requires an extra request. Result is cached."""
        if self.cid:
            results = get_json(self.cid, operation='sids')
            return results['InformationList']['Information'][0]['SID'] if results else []

    @memoized_property
    def aids(self):
        """Requires an extra request. Result is cached."""
        if self.cid:
            results = get_json(self.cid, operation='aids')
            return results['InformationList']['Information'][0]['AID'] if results else []

    @property
    def coordinate_type(self):
        if CoordinateType.TWO_D in self.record['coords'][0]['type']:
            return '2d'
        elif CoordinateType.THREE_D in self.record['coords'][0]['type']:
            return '3d'

    @property
    def charge(self):
        """Formal charge on this Compound."""
        return self.record['charge'] if 'charge' in self.record else 0

    @property
    def molecular_formula(self):
        """Molecular formula."""
        return _parse_prop({'label': 'Molecular Formula'}, self.record['props'])

    @property
    def molecular_weight(self):
        """Molecular Weight."""
        return _parse_prop({'label': 'Molecular Weight'}, self.record['props'])

    @property
    def canonical_smiles(self):
        """Canonical SMILES, with no stereochemistry information."""
        return _parse_prop({'label': 'SMILES', 'name': 'Canonical'}, self.record['props'])

    @property
    def isomeric_smiles(self):
        """Isomeric SMILES."""
        return _parse_prop({'label': 'SMILES', 'name': 'Isomeric'}, self.record['props'])

    @property
    def inchi(self):
        """InChI string."""
        return _parse_prop({'label': 'InChI', 'name': 'Standard'}, self.record['props'])

    @property
    def inchikey(self):
        """InChIKey."""
        return _parse_prop({'label': 'InChIKey', 'name': 'Standard'}, self.record['props'])

    @property
    def iupac_name(self):
        """Preferred IUPAC name."""
        # Note: Allowed, CAS-like Style, Preferred, Systematic, Traditional are available in full record
        return _parse_prop({'label': 'IUPAC Name', 'name': 'Preferred'}, self.record['props'])

    @property
    def xlogp(self):
        """XLogP."""
        return _parse_prop({'label': 'Log P'}, self.record['props'])

    @property
    def exact_mass(self):
        """Exact mass."""
        return _parse_prop({'label': 'Mass', 'name': 'Exact'}, self.record['props'])

    @property
    def monoisotopic_mass(self):
        """Monoisotopic mass."""
        return _parse_prop({'label': 'Weight', 'name': 'MonoIsotopic'}, self.record['props'])

    @property
    def tpsa(self):
        """Topological Polar Surface Area."""
        return _parse_prop({'implementation': 'E_TPSA'}, self.record['props'])

    @property
    def complexity(self):
        """Complexity."""
        return _parse_prop({'implementation': 'E_COMPLEXITY'}, self.record['props'])

    @property
    def h_bond_donor_count(self):
        """Hydrogen bond donor count."""
        return _parse_prop({'implementation': 'E_NHDONORS'}, self.record['props'])

    @property
    def h_bond_acceptor_count(self):
        """Hydrogen bond acceptor count."""
        return _parse_prop({'implementation': 'E_NHACCEPTORS'}, self.record['props'])

    @property
    def rotatable_bond_count(self):
        """Rotatable bond count."""
        return _parse_prop({'implementation': 'E_NROTBONDS'}, self.record['props'])

    @property
    def fingerprint(self):
        """Raw padded and hex-encoded fingerprint, as returned by the PUG REST API."""
        return _parse_prop({'implementation': 'E_SCREEN'}, self.record['props'])

    @property
    def cactvs_fingerprint(self):
        """PubChem CACTVS fingerprint.

        Each bit in the fingerprint represents the presence or absence of one of 881 chemical substructures.

        More information at ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
        """
        # Skip first 4 bytes (contain length of fingerprint) and last 7 bits (padding) then re-pad to 881 bits
        return '{0:020b}'.format(int(self.fingerprint[8:], 16))[:-7].zfill(881)

    @property
    def heavy_atom_count(self):
        """Heavy atom count."""
        if 'count' in self.record and 'heavy_atom' in self.record['count']:
            return self.record['count']['heavy_atom']

    @property
    def isotope_atom_count(self):
        """Isotope atom count."""
        if 'count' in self.record and 'isotope_atom' in self.record['count']:
            return self.record['count']['isotope_atom']

    @property
    def atom_stereo_count(self):
        """Atom stereocenter count."""
        if 'count' in self.record and 'atom_chiral' in self.record['count']:
            return self.record['count']['atom_chiral']

    @property
    def defined_atom_stereo_count(self):
        """Defined atom stereocenter count."""
        if 'count' in self.record and 'atom_chiral_def' in self.record['count']:
            return self.record['count']['atom_chiral_def']

    @property
    def undefined_atom_stereo_count(self):
        """Undefined atom stereocenter count."""
        if 'count' in self.record and 'atom_chiral_undef' in self.record['count']:
            return self.record['count']['atom_chiral_undef']

    @property
    def bond_stereo_count(self):
        """Bond stereocenter count."""
        if 'count' in self.record and 'bond_chiral' in self.record['count']:
            return self.record['count']['bond_chiral']

    @property
    def defined_bond_stereo_count(self):
        """Defined bond stereocenter count."""
        if 'count' in self.record and 'bond_chiral_def' in self.record['count']:
            return self.record['count']['bond_chiral_def']

    @property
    def undefined_bond_stereo_count(self):
        """Undefined bond stereocenter count."""
        if 'count' in self.record and 'bond_chiral_undef' in self.record['count']:
            return self.record['count']['bond_chiral_undef']

    @property
    def covalent_unit_count(self):
        """Covalently-bonded unit count."""
        if 'count' in self.record and 'covalent_unit' in self.record['count']:
            return self.record['count']['covalent_unit']

    @property
    def volume_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Shape', 'name': 'Volume'}, conf['data'])

    @property
    def multipoles_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Shape', 'name': 'Multipoles'}, conf['data'])

    @property
    def conformer_rmsd_3d(self):
        coords = self.record['coords'][0]
        if 'data' in coords:
            return _parse_prop({'label': 'Conformer', 'name': 'RMSD'}, coords['data'])

    @property
    def effective_rotor_count_3d(self):
        return _parse_prop({'label': 'Count', 'name': 'Effective Rotor'}, self.record['props'])

    @property
    def pharmacophore_features_3d(self):
        return _parse_prop({'label': 'Features', 'name': 'Pharmacophore'}, self.record['props'])

    @property
    def mmff94_partial_charges_3d(self):
        return _parse_prop({'label': 'Charge', 'name': 'MMFF94 Partial'}, self.record['props'])

    @property
    def mmff94_energy_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Energy', 'name': 'MMFF94 NoEstat'}, conf['data'])

    @property
    def conformer_id_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Conformer', 'name': 'ID'}, conf['data'])

    @property
    def shape_selfoverlap_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Shape', 'name': 'Self Overlap'}, conf['data'])

    @property
    def feature_selfoverlap_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Feature', 'name': 'Self Overlap'}, conf['data'])

    @property
    def shape_fingerprint_3d(self):
        conf = self.record['coords'][0]['conformers'][0]
        if 'data' in conf:
            return _parse_prop({'label': 'Fingerprint', 'name': 'Shape'}, conf['data'])



def get_compounds(identifier, namespace='cid', searchtype=None, as_dataframe=False, **kwargs):
    """Retrieve the specified compound records from PubChem.

    :param identifier: The compound identifier to use as a search query.
    :param namespace: (optional) The identifier type, one of cid, name, smiles, sdf, inchi, inchikey or formula.
    :param searchtype: (optional) The advanced search type, one of substructure, superstructure or similarity.
    :param as_dataframe: (optional) Automatically extract the :class:`~pubchempy.Compound` properties into a pandas
                         :class:`~pandas.DataFrame` and return that.
    """
    results = get_json(identifier, namespace, searchtype=searchtype, **kwargs)
    compounds = [Compound(r) for r in results['PC_Compounds']] if results else []
    if as_dataframe:
        return compounds_to_frame(compounds)
    return compounds

def compounds_to_frame(compounds, properties=None):
    """Construct a pandas :class:`~pandas.DataFrame` from a list of :class:`~pubchempy.Compound` objects.

    Optionally specify a list of the desired :class:`~pubchempy.Compound` properties.
    """
    import pandas as pd
    if isinstance(compounds, Compound):
        compounds = [compounds]
    properties = set(properties) | set(['cid']) if properties else None
    return pd.DataFrame.from_records([c.to_dict(properties) for c in compounds], index='cid')