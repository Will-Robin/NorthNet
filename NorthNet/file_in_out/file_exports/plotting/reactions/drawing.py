def draw_reaction(reaction_smiles,name):
    '''
    For drawing reactions
    Parameters
    ----------
    reaction_smiles: str
        Reaction SMILES for the reaction to be plotted.
    name: str
        Name for the output file.
    Returns
    -------
    fname: str
        Name of the file output.
    '''
    name = name.replace(".","_").replace(">>","__")
    DrawingOptions.atomLabelFontSize = 110
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 6.0

    reacs = [Chem.MolFromSmiles(x) for x in reaction_smiles.split(">>")[0].split(".")]
    prods = [Chem.MolFromSmiles(x) for x in reaction_smiles.split(">>")[1].split(".")]

    #reacs = [x for x in reacs if x.GetNumAtoms() > 1]
    #prods = [x for x in prods if x.GetNumAtoms() > 1]
    print(reacs)
    print(prods)
    for r in reacs:
        for atom in r.GetAtoms():
            atom.SetAtomMapNum(0)

    for p in prods:
        for atom in p.GetAtoms():
            atom.SetAtomMapNum(0)

    rxn = AllChem.ChemicalReaction() # Create a chemcical reaction
    [rxn.AddReactantTemplate(r) for r in reacs]
    [rxn.AddProductTemplate(p) for p in prods]

    #rxn = AllChem.ReactionFromSmarts(reaction_smiles)
    img = Draw.ReactionToImage(rxn)
    img.save('{}.png'.format(name))

    return '{}.png'.format(name)

def draw_reaction_SMARTS(reaction_smiles,name):
    '''
    For drawing reactions
    Parameters
    ----------
    reaction_smiles: str
        Reaction SMILES for the reaction to be plotted.
    name: str
        Name for the output file.
    Returns
    -------
    fname: str
        Name of the file output.
    '''
    name = name.replace(".","_").replace(">>","__")
    DrawingOptions.atomLabelFontSize = 110
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 6.0

    reacs = [Chem.MolFromSmarts(x) for x in reaction_smiles.split(">>")[0].split(".")]
    prods = [Chem.MolFromSmarts(x) for x in reaction_smiles.split(">>")[1].split(".")]

    #reacs = [x for x in reacs if x.GetNumAtoms() > 1]
    #prods = [x for x in prods if x.GetNumAtoms() > 1]

    for r in reacs:
        for atom in r.GetAtoms():
            atom.SetAtomMapNum(0)

    for p in prods:
        for atom in p.GetAtoms():
            atom.SetAtomMapNum(0)

    rxn = AllChem.ChemicalReaction() # Create a chemcical reaction
    [rxn.AddReactantTemplate(r) for r in reacs]
    [rxn.AddProductTemplate(p) for p in prods]

    #rxn = AllChem.ReactionFromSmarts(reaction_smiles)
    img = Draw.ReactionToImage(rxn)
    img.save('{}.png'.format(name))

    return '{}.png'.format(name)
