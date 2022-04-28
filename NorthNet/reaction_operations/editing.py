from NorthNet.text_parsing import conversions as conv
from NorthNet import molecule_operations as mol_ops


def canonicalise_reaction_smiles(reaction_smiles):
    """
    Canonicalise the consituent SMILES of a reaction SMILES string.

    Parameters
    ----------
    reaction_smiles: str

    Returns
    -------
    canonicalised_rxn_smiles: str
    """

    reactants, products = conv.reaction_smiles_split(reaction_smiles)

    canonical_reacts = [mol_ops.canonicalise(r) for r in reactants]
    canonical_prods = [mol_ops.canonicalise(p) for p in products]

    lhs = ".".join(canonical_reacts)
    rhs = ".".join(canonical_prods)

    canonicalised_rxn_smiles = f"{lhs}>>{rhs}"

    return canonicalised_rxn_smiles
