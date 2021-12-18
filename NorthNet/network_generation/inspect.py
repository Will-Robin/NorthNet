from itertools import compress

def reactive_species(species_list, substructures):
    '''
    TODO: rename

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
    # testing each molecule in the list of species for reaction group matches
    for mol in species_list: 
        for substructure in substructures:
            if mol.Mol.HasSubstructMatch(substructure):
                # append the molecule to the matches if
                # it contains the reactive group
                matches.append(mol) 
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
    for reaction in reactions:
        tag = True
        for exc in invalid_substructures:

            if any([p.HasSubstructMatch(exc.Mol)
                                            for p in reaction.Reaction.GetProducts()]):
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

def check_reaction_occurence(compound, network, reaction_template):
    '''
    Check if a reaction template has already been applied to a compound.
    
    Parameters
    ----------
    compound: NorthNet Compound object
    network: NorthNet Network object
    reaction_template: NorthNet ReactionTemplate object

    Returns
    -------
    bool
        Whether reaction type has been applied to the compound or not.
    '''

    reaction_name = reaction_template.Name

    reactions = []
    reactions.extend(compound.In)
    reactions.extend(compound.Out)

    used_reaction_classes = [network.get_reaction_name(r) for r in reactions]

    if reaction_name in used_reaction_classes:
        return True

    return False
