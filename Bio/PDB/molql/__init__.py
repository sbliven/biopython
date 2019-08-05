# Copyright 2019 by Spencer Bliven. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Manipulating and executing structural selection queries.

Allows writing queries in any of the MolQL-supported query syntaxes.

See http://molql.org/
"""
import warnings
import re

try:
    import pyparsing as pp
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError("Install pyparsing to use MolQL "
                                       "(e.g. pip install pyparsing)")

from .. import BiopythonWarning
from . import Select
from . import extract_to_structure


_hydrogen = re.compile("[123 ]*H.*")


# "range" format
class ResRange(object):
    """Represents a simple residue range

    Ranges are inclusive between start and end. If one endpoint is missing
    than the beginning or end of the chain is implied

    Arguments:
    - chain_id (str): Author chain ID
    - start (int or (int,str)): starting residue number, optionally with insertion code
    - end (int or (int,str)): ending residue number, optionally with insertion code
    """
    def __init__(self, chain_id, start=None, end=None):
        self.chain_id = chain_id
        if type(start) is tuple:
            self.start = start
        else:
            self.start = (start, " ")
        if type(end) is tuple:
            self.end = end
        else:
            self.end = (end, " ")

    def __str__(self):
        return "%s.%s%s-%s%s" % (self.chain_id,
                                 "^" if self.start is None or self.start[0] is None else self.start[0],
                                 "" if self.start is None or self.start[1] is None else self.start[1].strip(),
                                 "$" if self.end is None or self.end[0] is None else self.end[0],
                                 "" if self.end is None or self.end[1] is None else self.end[1].strip())

    def __repr__(self):
        return str(self)


_ranges_bnf = None


def rangesBNF():
    """Construct the grammar for the "range" format

    Returns a grammar, which upon successfully parsing should produce a list of ResRange.

    """
    global _ranges_bnf
    if _ranges_bnf is None:
        def parseResi(s, l, t):
            if len(t) == 1:
                return (t[0], " ")
            else:
                return tuple(t)

        def parseRange(s, l, t):
            return ResRange(*t)

        inscode = pp.Word(pp.alphanums)
        signedInt = pp.Regex("[+-]?[0-9]+").setParseAction(lambda s, l, t: [int(t[0])])
        resi = (signedInt + pp.Optional(inscode)).leaveWhitespace().setParseAction(parseResi)
        chain_id = pp.Word(pp.alphanums)
        rangeExpr = (chain_id.setResultsName('chain_id') +
                     pp.Optional(pp.Suppress('.') + resi.setResultsName('start') +
                                 pp.Optional(pp.Suppress('-') + resi.setResultsName('end')))
                     ).leaveWhitespace().setParseAction(parseRange)
        ranges = pp.delimitedList(rangeExpr)
        _ranges_bnf = ranges
    return _ranges_bnf


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
            self.resRanges = rangesBNF().parseString(query)
        else:
            raise ValueError("Unsupported query format: %s" % format)

    def __call__(self, structure):
        """Evaluate query on a structure

        Args:
        - structure (Bio.PDB.Structure): structure to filter

        Returns:
        Bio.PDB.Structure: result of the query
        """
        if self.resRanges is None:
            return structure
        else:

            return extract_to_structure(structure, ChainSelector(self.resRanges))
