def reconstitute_reaction(reactants,products):

    '''
    Parameters
    ----------
    reactants: list
        List of reactants.
    products: list
        List of products
    Returns
    -------
    reaction_string: str
        Reaction SMILES
    '''

    LHS = ".".join(reactants)
    RHS = ".".join(products)

    reaction_string = LHS + ">>" + RHS

    return reaction_string

def species_from_reaction(reactions):
    '''Take a set of reactions and return prducts and reactants
    as a list of species.'''
    species = []
    for u in reactions:
        us = u.split('>>')
        for s in us:
            sp = s.split('.')
            species.extend(sp)

    species = list(set(species))

    return species

def get_product_smiles(reaction):
    '''
    Input
    reaction: NorthNet Generated_Reaction
    Output
    smiles: str
    '''
    prods = reaction.Reaction.GetProducts()
    prods = [remove_atom_mapping(p) for p in prods]
    smiles = [Chem.MolToSmiles(p) for p in prods]
    return smiles
