# Copyright 2019 by Spencer Bliven. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.MolQL package"""
import unittest
import json
from io import BytesIO
from Bio.PDB.molql import MolQL
from Bio.PDB.molql.range_script import rangesBNF, ResRange
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Dice


class MolQLTests(unittest.TestCase):
    """Tests for PDB.MolQL module."""

    def test_ranges_parsing(self):
        """Test parsing of ranges directly."""
        import pyparsing as pp

        tree = rangesBNF().parseString('A', True)
        self.assertEqual(len(tree), 1)
        self.assertEqual(tree[0].chain_id, 'A')
        self.assertEqual(tree[0].start, (None, " "))
        self.assertEqual(tree[0].end, (None, " "))

        tree = rangesBNF().parseString('A.1-10,B,C.1', True)
        self.assertEqual(len(tree), 3)
        self.assertEqual(tree[0].chain_id, "A")
        self.assertEqual(tree[0].start, (1, " "))
        self.assertEqual(tree[0].end, (10, " "))
        self.assertEqual(tree[1].chain_id, 'B')
        self.assertEqual(tree[1].start, (None, " "))
        self.assertEqual(tree[1].end, (None, " "))
        self.assertEqual(tree[2].chain_id, 'C')
        self.assertEqual(tree[2].start, (1, " "))
        self.assertEqual(tree[2].end, (None, " "))

        with self.assertRaises(pp.ParseException):
            rangesBNF().parseString('1-2-3', True)

    def test_ranges_dicing(self):
        """Test MolQL dicing from ranges format

        """

        # Use the same test cases as test_PDB_Dice
        pdb_file = "PDB/2BEG.pdb"
        parser = PDBParser()
        structure = parser.get_structure("scr", pdb_file)

        # Substructure
        query = MolQL("B.18-20", "range")
        substructure = query(structure)

        self.assertEqual(sum(1 for c in substructure.get_chains()), 1)
        self.assertEqual(sum(1 for r in substructure.get_residues()), 3)
        self.assertEqual(sum(1 for a in substructure.get_atoms()), 29)

        # Write reference using extract for comparison
        file_pdb_extract = BytesIO()
        Dice.extract(structure, "B", 18, 20, file_pdb_extract)

        # Test it matches the MolQL version
        io = PDBIO()
        io.set_structure(substructure)
        file_pdb_molql = BytesIO()
        io.save(file_pdb_molql)

        # Check files equal
        self.assertEqual(len(file_pdb_molql.getvalue()),
            len(file_pdb_extract.getvalue()))
        self.assertEqual(file_pdb_molql.getvalue(), file_pdb_extract.getvalue())

    def test_json(self):

        json_example2 = """{
            "source": "molql-explorer",
            "version": "0.1.0",
            "expression": {
                "head": "structure.generator.atom-groups",
                "args": {
                    "residue-test": {
                        "head": "core.rel.eq",
                        "args": {
                            "0": {
                                "head": "structure.atom-property.macromolecular.label_comp_id"
                            },
                            "1": "ALA"
                        }
                    }
                }
            }
        }"""

        query = MolQL(json_example2, "json")

        # Test writing
        j = query.dumps("json")

        expected_dict = json.loads(json_example2)
        result_dict = json.loads(j)

        # ignore source
        del expected_dict["source"]
        del result_dict["source"]

        self.assertDictEqual(result_dict, expected_dict)

        json_example1 = """{
            "source": "molql-explorer",
            "version": "0.1.0",
            "expression": {
                "head": "structure.generator.atom-groups",
                "args": {
                    "residue-test": {
                        "head": "core.rel.eq",
                        "args": {
                            "0": {
                                "head": "structure.atom-property.macromolecular.auth_comp_id"
                            },
                            "1": "ALA"
                        }
                    },
                    "atom-test": {
                        "head": "core.set.has",
                        "args": {
                            "0": {
                                "head": "core.type.set",
                                "args": {
                                    "0": {
                                        "head": "structure.type.element-symbol",
                                        "args": {
                                            "0": "C"
                                        }
                                    },
                                    "1": {
                                        "head": "structure.type.element-symbol",
                                        "args": {
                                            "0": "N"
                                        }
                                    }
                                }
                            },
                            "1": {
                                "head": "structure.atom-property.core.element-symbol"
                            }
                        }
                    }
                }
            }
        }
        """

        # Test parsing
        query = MolQL(json_example1, "json")

        # Test writing
        j = query.dumps("json")

        expected_dict = json.loads(json_example1)
        result_dict = json.loads(j)

        # ignore source
        del expected_dict["source"]
        del result_dict["source"]

        self.assertDictEqual(result_dict, expected_dict)


    def test_molql_lisp(self):
        lisp_example1 = """(structure.generator.atom-groups
  :residue-test (core.rel.eq
    (structure.atom-property.macromolecular.auth_comp_id)
    ALA)
  :atom-test (core.set.has
    (core.type.set
      (structure.type.element-symbol C)
      (structure.type.element-symbol N))
    (structure.atom-property.core.element-symbol)))
"""


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
