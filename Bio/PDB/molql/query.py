"""Execute a MolQL query on BioPython objects
"""
import warnings
import re
import abc
from functools import wraps

from ... import BiopythonWarning
from .. import Select
from .. import extract_to_structure
from . import syntax
# Make sure all known syntax modules are loaded
from . import json_script, molql_script, range_script
# Load known symbols
from . import core, structure
from . import symbols


class ChainSelector(Select):
    """
    Only accepts residues with right chainid
    and between start and end. Remove hydrogens, waters and ligands.
    Only use model 0 by default.
    """
    _hydrogen = re.compile("[123 ]*H.*")

    def __init__(self, ranges):
        """"""
        self.ranges = ranges

    def accept_chain(self, chain):
        if any(chain.get_id() == r.chain_id for r in self.ranges):
            return 1
        return 0

    def accept_residue(self, residue):
        # residue - between start and end
        hetatm_flag, resseq, icode = residue.get_id()
        chain_id = residue.get_parent().get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode is not None and icode.strip() != "":
            warnings.warn("WARNING: Icode %s at position %s" %
                          (icode, resseq), BiopythonWarning)
        if any(chain_id == r.chain_id and r.start[0] <= resseq <= r.end[0]
               for r in self.ranges):
            return 1
        return 0

    def accept_atom(self, atom):
        """Verify if atoms are not Hydrogen."""
        # atoms - get rid of hydrogens
        name = atom.get_id()
        if ChainSelector._hydrogen.match(name):
            return 0
        else:
            return 1


class StructureContext(object):
    """Tracks the application of a MolQL query to a single structure
    """

    def __init__(self, query, structure):
        self.query = query
        self.structure = structure

    def _eval(self, expr):
        """Evaluates an Expression in this context

        Looks up implementations in the symbol table and evaluate them recursively.
        """
        if not isinstance(expr, syntax.Expression):
            return expr

        impl = symbols.get_symbol(expr.head)
        posargs = {}
        kwargs = {}
        for k, v in expr.args:
            try:
                posargs[int(k)] = v
            except ValueError:
                kwargs[k] = v

        if len(posargs) > 0:
            maxarg = max(posargs.keys())
            # Raises KeyError for non-sequential positional arguments
            args_evaled = (self._eval(posargs[i]) for i in xrange(maxarg))
        else:
            args_evaled = []
        kwargs_evaled = {k: self._eval(v) for k, v in kwargs}
        return impl(self, *args_evaled, **kwargs_evaled)

    def __call__(self):
        return self._eval(self.query.expression)


class MolQL(object):
    """Represent a molecular query. Queries can be created from a variety
    of selection syntaxes. Queries can then be applied to Bio.PDB.Structure
    objects to reduce to the selected substructure.

    In their most general form, MolQL queries act on a set of structure
    fragments (conceptually a list of Atom lists).

    Supported formats
    -----------------
    MolQL can be extended with additional formats. The complete list of formats
    can be accessed with `molql.syntax.get_syntaxes()`. Some standard formats:

    - *auto* Auto-detect the format (not recommended)
    - *range* Simple residue range format. Example: "A.1-10"
    - *molql* MolQL Lisp-like
    - *json* MolQL JSON
    - *pymol* Pymol syntax
    - *jmol* Jmol syntax
    - *vmd* VMD syntax

    Examples
    --------
    >>> from Bio.PDB import MolQL
    >>> from Bio.PDB import PDBParser
    >>> parser = PDBParser()
    >>> structure = parser.get_structure("2BEG", "PDB/2BEG.pdb")
    >>> query = MolQL("B.18-20", "range")
    >>> substructure = query(structure)
    >>> list(substructure.get_residues())
    [<Residue VAL het=  resseq=18 icode= >,
     <Residue PHE het=  resseq=19 icode= >,
     <Residue PHE het=  resseq=20 icode= >]

    """

    def __init__(self, query, format="auto"):
        """Create a new query object.

        Typical use is to specify the query and format. In this case the query
        is parsed according to the
        Args:
        - query (str): selection query
        - format (str): Format of the query (case insensitive)
        """
        try:
            parser = syntax.get_syntax(format)
        except KeyError:
            raise ValueError("Unsupported query format: %s" % format)

        self.expression = parser.loads(query)

    def dice(self, structure):
        """Evaluate query on a structure

        Returns a set of Structures representing each fragment of the result

        Args:
        - structure (Bio.PDB.Structure): structure to filter

        Returns:
        Bio.PDB.Structure: result of the query
        """
        raise NotImplementedError()

    def apply(self, structure):
        """Low-level function to evaluate a query on a set of fragments.

        Args:
        - structure (Bio.PDB.Structure): structure to filter

        Returns:
        Iterable[Iterable[Atom]]: Set of all fragments
        """
        context = StructureContext(self, structure)
        return context()

    def __call__(self, structure):
        """Evaluate query on a structure

        If multiple fragments would be returned by the query, performs an
        implicit union over all fragments to return a single structure.

        Args:
        - structure (Bio.PDB.Structure): structure to filter

        Returns:
        Bio.PDB.Structure: result of the query
        """
        # TODO stub
        if self.expression is None:
            return structure
        else:
            # TODO Broken! 'ranges' misuses the expression to store a list of ResRange
            return extract_to_structure(structure, ChainSelector(self.expression))

    def dump(self, fp, format, **kw):
        try:
            parser = syntax.get_syntax(format)
        except KeyError:
            raise ValueError("Unsupported query format: %s" % format)

        parser.dump(self.expression, fp, **kw)

    def dumps(self, format, **kw):
        try:
            parser = syntax.get_syntax(format)
        except KeyError:
            raise ValueError("Unsupported query format: %s" % format)

        return parser.dumps(self.expression, **kw)

    def __str__(self):
        return self.dumps('molql')