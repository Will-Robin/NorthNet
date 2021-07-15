from itertools import compress

def reactive_species(species_list, substructures):
    '''
    Find species with the given subtructure.
    Returns a list of matching compounds

    Parameters
    ----------
    species_list: list
        List of NorthNet Compound objects from which
        reactive molecules are extracted.

    substructure: list of rdkit mol objects
        Reactive substructures.

    Returns
    -------
    matches: list
        List of NorthNet Compound objects which contain the substructure.
    '''
    matches = [] # empty container to put matching molecules in
    for mol in species_list: # testing each molecule in the list of species for reaction group matches
        for s in substructures:
            if mol.Mol.HasSubstructMatch(s):
                matches.append(mol) # append the molecule to the matches if in contains the reactive group
    return matches

def remove_invalid_reactions(reactions,invalid_substructures):

    '''
    Designed to remove reactions with products which contain
    invalid substructures as defined by the user.

    Parameters
    ----------
    reactions: list
        List of NorthNet Generated_Reaction objects.
    invalid_substructures: list
        list of NorthNet Substructure objects.

    Returns
    -------
    reactions: list
        list of reactions with those that produce
        invalid substructures removed.
    '''

    sortlist = []
    for r in reactions:
        tag = True
        for exc in invalid_substructures:
            if any([p.HasSubstructMatch(exc.Mol)
                                            for p in r.Reaction.GetProducts()]):
                tag = False

        sortlist.append(tag)

    # new list with compounds tagged false removed
    reactions = list(compress(reactions, sortlist))

    return reactions

def check_reaction_input(reactant_list, reactive_substructs):
    '''
    Checks for valid reaction input by checking the reactants
    contain the supplied substructures and are in the same order.

    Returns False if reaction inputs are incompatible with
    the reactive substructures provided.

    reactant_list: NorthNet Compound Object

    reaction_template: NorthNet ReactionTemplate object
    '''

    test_list = [r.Mol.HasSubstructMatch(s)
                            for r,s in zip(reactant_list,reactive_substructs)]

    return all(test_list)
