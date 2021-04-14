def balance_reaction(reaction_object, element):
    '''
    Balances a reaction by counting the number of element atoms on either side
    of it and adding elements on the other side if the stoichiometry is unbalanced.

    Parameters
    ----------
    reaction_object: NorthNet Reaction
        Reaction Object to be balanced.
    element: int or str
        Atomic number (int) or SMARTS for element to be balanced.
    Returns
    -------
    Modifies reaction_object
    '''

    res = reaction_object.Reactants
    prs = reaction_object.Products

    if reaction_object.ReactionSMILES.count(element) == 0:
        pass
    else:
        if res.count(element) != prs.count(element):
            if res.count(element) > prs.count(element):
                n = res.count(element) - prs.count(element)
                prs.extend([element]*n)

            if res.count(element) < prs.count(element):
                n = prs.count(element) - res.count(element)
                res.extend([element]*n)

            reaction_object.update_reaction()

        else:
            pass

def cleanup(reactions):
    '''Used to clear up the output of rdkit reaction function, parsing multiple
    reaction outcomes. Not perfect.'''

    reactions = [r.replace('[H+][O-]','O')    for r in reactions]
    reactions = [r.replace('[O-][H+]','O')    for r in reactions]
    #reactions = [r.replace('=[C@H]'   ,'=C')  for r in reactions]
    #reactions = [r.replace('=[C@@H]'   ,'=C') for r in reactions]
    reactions = [r.replace("/"   ,"") for r in reactions]
    reactions = [r.replace( r"\\"  ,"") for r in reactions]
    #reactions = [r.replace('@'   ,'') for r in reactions]

    reactions_out = []
    for c in range(0,len(reactions)):
        # Split the reaction up for inspection
        LHS = reactions[c].split(">>")[0].split(".")
        RHS = reactions[c].split(">>")[1].split(".")
        ins = ''
        for l in LHS:
            mol = Chem.MolFromSmiles(l)
            x = Chem.MolToSmiles(mol, canonical = True)
            ins += x + '.'

        ins = ins.strip('.')
        ins += '>>'

        for r in RHS:
            mol = Chem.MolFromSmiles(r)
            x = Chem.MolToSmiles(mol, canonical = True)
            ins += x + '.'

        reactions_out.append(ins.strip('.'))

    return reactions_out

def remove_reaction_mapping(reaction):
    '''
    Input
    reaction: NorthNet Generated_Reaction
    Output
    mod_rxn: NorthNet Generated_Reaction
    '''
    reacs = reaction.Reaction.GetReactants()
    prods = reaction.Reaction.GetProducts()
    for r in reacs:
        for atom in r.GetAtoms():
            atom.SetAtomMapNum(0)

    for p in prods:
        for atom in p.GetAtoms():
            atom.SetAtomMapNum(0)

    rxn = AllChem.ChemicalReaction() # Create a chemcical reaction
    [rxn.AddReactantTemplate(r) for r in reacs]
    [rxn.AddProductTemplate(p) for p in prods]

    mod_rxn = Classes.Generated_Reaction(rxn, reaction.Reaction_Template)

    return mod_rxn
