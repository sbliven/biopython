# Copyright 2019 by Sergio Valqui. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.Dice Module."""

import unittest
import tempfile
import os
from io import BytesIO

from Bio.PDB import Dice
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO


class DiceTests(unittest.TestCase):
    """Tests for PDB.Dice module."""

    def test_dice(self):
        """Self test for PDB.Dice module."""
        # Load and dice 2BEG.pdb
        parser = PDBParser()
        pdb_file = "PDB/2BEG.pdb"
        structure = parser.get_structure("scr", pdb_file)
        file_number, file_pdb_diced = tempfile.mkstemp()
        os.close(file_number)
        # dice 2BEG.pdb, chain "B" from residue 18 to 20"
        Dice.extract(structure, "B", 18, 20, file_pdb_diced)
        try:
            # Load the diced file
            diced_str = parser.get_structure("2beg_diced", file_pdb_diced)
            l_diced_cha = list(diced_str.get_chains())
            l_diced_res = list(diced_str.get_residues())
            l_diced_ato = list(diced_str.get_atoms())

            # Chain checks
            self.assertEqual(len(l_diced_cha), 1)
            self.assertEqual(l_diced_cha[0].id, "B")

            # Residue Checks
            self.assertEqual(len(l_diced_res), 3)
            self.assertEqual(l_diced_res[0].id[1], 18)
            self.assertEqual(l_diced_res[2].id[1], 20)

            # Atom checks
            self.assertEqual(len(l_diced_ato), 29)
            self.assertEqual(l_diced_ato[0].name, 'N')
            self.assertEqual(l_diced_ato[0].parent.resname, 'VAL')
            self.assertEqual(l_diced_ato[0].parent.parent.id, 'B')
            self.assertEqual(l_diced_ato[28].name, 'CZ')
            self.assertEqual(l_diced_ato[28].parent.resname, 'PHE')
            self.assertEqual(l_diced_ato[28].parent.parent.id, 'B')

        finally:
            os.remove(file_pdb_diced)

    def test_extract_to_structure(self):
        "Test substructure extraction"

        pdb_file = "PDB/2BEG.pdb"
        parser = PDBParser()
        structure = parser.get_structure("scr", pdb_file)

        # Substructure
        sel = Dice.ChainSelector("B", 18, 20)
        substructure = Dice.extract_to_structure(structure, sel)

        self.assertEqual(sum(1 for c in substructure.get_chains()), 1)
        self.assertEqual(sum(1 for r in substructure.get_residues()), 3)
        self.assertEqual(sum(1 for a in substructure.get_atoms()), 29)

        # Write reference using extract
        file_pdb_extract = BytesIO()
        Dice.extract(structure, "B", 18, 20, file_pdb_extract)

        # Test it matches the extract_from_structure
        io = PDBIO()
        io.set_structure(substructure)
        file_pdb_extract_from_struct = BytesIO()
        io.save(file_pdb_extract_from_struct)

        # Check files
        # with open("extract.pdb", "wb") as f1:
        #     f1.write(file_pdb_extract.getvalue())
        # with open("extract_struct.pdb", "wb") as f2:
        #     f2.write(file_pdb_extract_from_struct.getvalue())
        self.assertEqual(len(file_pdb_extract.getvalue()),
            len(file_pdb_extract_from_struct.getvalue()))
        self.assertEqual(file_pdb_extract.getvalue(), file_pdb_extract_from_struct.getvalue())


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
