# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for chopping up (dicing) a structure.

This module is used internally by the Bio.PDB.extract() function.
"""

import re
import warnings

from Bio.PDB.PDBIO import PDBIO, Select
from Bio import BiopythonWarning
from Bio.PDB.StructureBuilder import StructureBuilder

_hydrogen = re.compile("[123 ]*H.*")


class ChainSelector(Select):
    """Only accepts residues with right chainid, between start and end.

    Remove hydrogens, waters and ligands. Only use model 0 by default.
    """

    def __init__(self, chain_id, start, end, model_id=0):
        """Initialize the class."""
        self.chain_id = chain_id
        self.start = start
        self.end = end
        self.model_id = model_id

    def accept_model(self, model):
        """Verify if model match the model identifier."""
        # model - only keep model 0
        if model.get_id() == self.model_id:
            return 1
        return 0

    def accept_chain(self, chain):
        """Verify if chain match chain identifier."""
        if chain.get_id() == self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # residue - between start and end
        hetatm_flag, resseq, icode = residue.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            warnings.warn("WARNING: Icode %s at position %s"
                          % (icode, resseq), BiopythonWarning)
        if self.start <= resseq <= self.end:
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


def extract(structure, chain_id, start, end, filename):
    """Write out selected portion to filename."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)


def extract_to_structure(structure, select=Select()):
    """Construct a new Structure that is a subset of the given structure

    Args:
    - structure (Bio.PDB.Structure): input structure
    - select (Bio.PDB.PDBIO.Select): Filter instance specifying whether to
      include each model, chain, residue, and atom
    """
    builder = StructureBuilder()
    current_segid = None

    def build_atom(atom, builder=builder, current_segid=current_segid):
        if select.accept_atom(atom):
            builder.init_atom(
                atom.get_name(),
                atom.get_coord(),
                atom.get_bfactor(),
                atom.get_occupancy(),
                atom.get_altloc(),
                atom.get_fullname(),
                atom.get_serial_number(),
                atom.element
            )

    def build_residue(res, builder=builder, current_segid=current_segid):
        if select.accept_residue(res):
            segid = res.get_segid()
            if segid != current_segid:
                current_segid = segid
                builder.init_seg(current_segid)
            builder.init_residue(res.get_resname(), *res.get_id())
            for atom in res:
                build_atom(atom)

    def build_chain(chain, builder=builder, current_segid=current_segid):
        if select.accept_chain(chain):
            builder.init_chain(chain.get_id())
            for res in chain:
                build_residue(res)

    def build_model(model, builder=builder, current_segid=current_segid):
        if select.accept_model(model):
            builder.init_model(model.get_id(), model.serial_num)
            for chain in model:
                build_chain(chain)

    builder.init_structure(structure.get_id())

    for model in structure:
        build_model(model)

    builder.set_header(structure.header)

    return builder.get_structure()


if __name__ == "__main__":

    from Bio.PDB.PDBParser import PDBParser

    import sys

    p = PDBParser()
    s = p.get_structure("scr", sys.argv[1])

    extract(s, " ", 1, 100, "out.pdb")
