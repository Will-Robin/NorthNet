from rdkit import Chem
from rdkit.Chem import AllChem


def canonicalise(smiles):
    """
    Parameters
    ----------
    smiles: str
        Canonicalises a SMILES string.

    Return
    ------
    str
        Canonicalised smiles.

    """
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        print(f"Did not canonicalise {smiles}.")
        return smiles


def mirror_smiles(smiles):
    """
    Reverses enantiomers via string replacement in SMILES

    Parameters
    ----------
    smiles: str
        SMILES string to be mirrored.
    Returns
    -------
    m1: str
        Mirrored smiles
    """

    molecule_1 = smiles.replace("@@", "ZY")
    molecule_2 = molecule_1.replace("@", "@@")
    molecule_3 = molecule_2.replace("ZY", "@")
    canonicalised_molecule = canonicalise(molecule_3)

    return canonicalised_molecule


def mirror(compound):
    """
    Reverses the chirality of a NorthNet Compound object.

    Parameters
    ----------
    compound: NorthNet.Classes.Compound
        Compound to be mirrored. Modifed in place.
    Returns
    -------
    None
    """

    m1 = mirror_smiles(compound.SMILES)

    compound.Mol = Chem.MolFromSmiles(m1)
    compound.SMILES = Chem.MolToSmiles(compound.Mol, isomericSmiles=True)


def incorrect_chiral_H_solve(mol):
    """
    Checks for carbon centres which exceed a total valence of 4, sets
    their chirality to unspecified and removes explicit hydrogens.

    Example input:
    OC[C@H](O)C(O)=[C@H](O)[C@H](O)CO
                    ^- This carbon should not contain chiral information.

    Example output given input above:
    OC[C@H](O)C(O)=C(O)[C@H](O)CO
                   ^- Chiral information has been removed.

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol
        Molecule to be corrected.
    Returns
    -------
    mol: rdkit.Chem.rdchem.Mol
        Cleaned molecule
    """
    mol.UpdatePropertyCache(strict=False)
    Chem.AssignStereochemistry(
        mol, cleanIt=True, force=True, flagPossibleStereoCenters=True
    )

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetTotalValence() > 4:
            atom.SetChiralTag(AllChem.CHI_UNSPECIFIED)
            atom.SetNumExplicitHs(0)

    return mol
