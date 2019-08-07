"""Execute a MolQL query on BioPython objects
"""
import warnings
import re

from ... import BiopythonWarning
from .. import Select
from .. import extract_to_structure
from . import ranges_script

_hydrogen = re.compile("[123 ]*H.*")


class ChainSelector(Select):
    """
    Only accepts residues with right chainid
    and between start and end. Remove hydrogens, waters and ligands.
    Only use model 0 by default.
    """
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
        if _hydrogen.match(name):
            return 0
        else:
            return 1


class MolQL(object):
    """Represent a molecular query. Queries can be created from a variety
    of selection syntaxes. Queries can then be applied to Bio.PDB.Structure
    objects to reduce to the selected substructure.

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

        Supported formats:
        - *range* Simple residue range format. Example: "A.1-10"

        Planned formats:
        - *auto* Auto-detect format (not recommended)
        - *molql* MolQL Lisp-like
        - *json* MolQL JSON
        - *pymol* Pymol syntax
        - *jmol* Jmol syntax
        - *vmd* VMD syntax

        Args:
        - query (str): selection query
        - format (str): Format of the query (case insensitive)
        """
        if format.lower() == "range":
            # TODO convert to expression
            self.resRanges = ranges_script.rangesBNF().parseString(query)
        elif format.lower() == "json":
            self.expression = None  # TODO
        else:
            raise ValueError("Unsupported query format: %s" % format)

    def __call__(self, structure):
        """Evaluate query on a structure

        Args:
        - structure (Bio.PDB.Structure): structure to filter

        Returns:
        Bio.PDB.Structure: result of the query
        """
        # TODO stub
        if self.resRanges is None:
            return structure
        else:

            return extract_to_structure(structure, ChainSelector(self.resRanges))
