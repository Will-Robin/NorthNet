from rdkit import Chem
from rdkit.Chem import AllChem

def canonicalise(smiles):
    '''
    Canonicalises a SMILES string
    '''
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles = True)

def mirror_smiles(smiles):
    '''
    Reverses enantiomers via string replacement in SMILES

    Parameters
    ----------
    smiles: str
        SMILES string to be mirrored.
    Returns
    -------
    m1: str
        Mirrored smiles
    '''

    m1 = smiles.replace("@@","ZY")
    m1 = m1.replace("@","@@")
    m1 = m1.replace("ZY","@")
    m1 = canonicalise(m1)

    return m1

def mirror(compound):
    '''
    Reverses the chirality of a NorthNet Compound object.

    Parameters
    ----------
    compound: NorthNet compound object
        Compound to be mirrored. Modifed in place.
    Returns
    -------
    None
    '''

    m1 = mirror_smiles(compound.SMILES)

    compound.Mol = Chem.MolFromSmiles(m1)
    compound.SMILES = Chem.MolToSmiles(compound.Mol, isomericSmiles = True)

def incorrect_chiral_H_solve(mol):
    '''
    Checks for carbon centres which exceed a total valence of 4, sets
    their chirality to unspecified and removes explicit hydrogens.

    Parameters
    ----------
    mol: rdkit mol object
        Molecule to be corrected.
    Returns
    -------
    mol: rdkit mol object
        Cleaned molecule
    '''
    mol.UpdatePropertyCache(strict = False)
    Chem.AssignStereochemistry(mol,
                            cleanIt=True,
                            force=True,
                            flagPossibleStereoCenters=True
                            )

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetTotalValence() > 4:
            atom.SetChiralTag(AllChem.CHI_UNSPECIFIED)
            atom.SetNumExplicitHs(0)

    return mol
