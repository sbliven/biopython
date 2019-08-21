"""Implements MolQL queries relating to protein structure on BioPython objects

"""
import itertools

from .symbols import symbol


@symbol("structure.generator.atom-groups")
def atom_groups(context, fragments=None, model_test=None, chain_test=None,
                residue_test=None, atom_test=None):
    """
    Args:
    - model_test (Model -> bool): Default: accept model 1
    - chain_test (Chain -> bool)
    - residue_test (Residue -> bool)
    - atom_test (Atom -> bool)

    Returns:
    Generator yielding Iterator[Atom]. Each atom is returned in its own list.
    """
    structure = context.structure
    # ignore input fragments and iterate over the whole structure
    models = structure.get_list()
    if model_test is None:
        models = itertools.islice(models, 1)
    for model in structure:
        if model_test is not None and not model_test(context, model):
            continue
        for chain in model:
            if chain_test is not None and not chain_test(context, chain):
                continue
            for residue in chain:
                # TODO deal with disordered residues?
                if residue_test is not None and not residue_test(context, residue):
                    continue
                for atom in residue:
                    if atom_test is not None and not atom_test(context, atom):
                        continue
                    yield [atom]
