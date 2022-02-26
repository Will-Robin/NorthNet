from rdkit import Chem
from rdkit.Chem import AllChem


def reaction_smiles_split(reaction_smiles):
    """
    Parameters
    ----------
    reaction_smiles: str


    Returns
    -------
    reactants, products: tuple(list[str])
    """

    split_rxn = reaction_smiles.split(">>")

    reactants = split_rxn[0].split(".")
    products = split_rxn[1].split(".")

    return reactants, products


def smiles_to_rdkit_reaction(smiles):
    """
    For converting reaction SMILES strings into RDKit AllChem.ChemicalReaction
    objects.

    Parameters
    ----------
    smiles: str
        Reaction SMILES.

    Returns
    -------

    rdkit_reaction: RDKit AllChem.Reaction
        Converted reaction.
    """

    products_reactants = smiles.split(">>")

    reactants = products_reactants[0].split(".")
    products = products_reactants[1].split(".")

    reactants_as_mol = [Chem.MolFromSmiles(r) for r in reactants]
    products_as_mol = [Chem.MolFromSmiles(p) for p in products]

    rdkit_reaction = AllChem.ChemicalReaction()

    [rdkit_reaction.AddReactantTemplate(r) for r in reactants_as_mol]
    [rdkit_reaction.AddProductTemplate(p) for p in products_as_mol]

    return rdkit_reaction
