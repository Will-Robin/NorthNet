def canonicalise(smiles):
    '''
    Canonicalises a SMILES molecular structure.
    '''
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles = True)

def mirror(compound):
    '''
    Reverses enantiomers

    Parameters
    ----------
    compound: NorthNet/ NetGen compound object
        Compound to be mirrored.
    Returns
    -------
    None
    '''

    m1 = compound.SMILES.replace("@@","ZY")
    m1 = m1.replace("@","@@")
    m1 = m1.replace("ZY","@")

    compound.Mol = Chem.MolFromSmiles(m1)
    compound.SMILES = Chem.MolToSmiles(compound.Mol, isomericSmiles = True)

def mirror_smiles(smiles):
    '''
    Reverses enantiomers

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

    return m1

def incorrect_chiral_H_solve(mol):
    '''
    Removes explicit hydrogens left over from deprotonation of a stereogenic carbon.
    Checks for carbon centres which exceed a total valence of 4 and those with
    explicit
    Parameters
    ----------
    mol: rdkit mol object
        Molecule to be corrected.
    Returns
    -------
    None
    '''
    mol.UpdatePropertyCache(strict = False)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True,flagPossibleStereoCenters=True)

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetTotalValence() > 4:
            atom.SetChiralTag(AllChem.CHI_UNSPECIFIED)
            atom.SetNumExplicitHs(0)

    return mol
