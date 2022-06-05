from itertools import compress


def get_reactive_compounds(species_list, substructures):
    """
    Find species with the given subtructure.
    Returns a list of matching compounds

    Parameters
    ----------
    species_list: list[NorthNet.Classes.Compound]
        List of NorthNet Compound objects from which
        reactive molecules are extracted.
    substructure: list[rdkit.Chem.rdchem.Mol]
        Reactive substructures.

    Returns
    -------
    matches: list[NorthNet.Classes.Compound]
        List of NorthNet Compound objects which contain the substructure.
    """
    matches = []
    # testing each molecule in the list of species for reaction group matches
    for mol in species_list:
        for substructure in substructures:
            if mol.Mol.HasSubstructMatch(substructure):
                matches.append(mol)
    return matches


def remove_invalid_reactions(reactions, invalid_substructures):
    """
    Removes reactions with products that contain invalid substructures.

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
    """

    sortlist = []
    for reaction in reactions:
        tag = True
        for exc in invalid_substructures:

            products = reaction.Reaction.GetProducts()

            substruct_matches = [p.HasSubstructMatch(exc.Mol) for p in products]

            if any(substruct_matches):
                tag = False

        sortlist.append(tag)

    # new list with compounds tagged as false removed
    reactions = list(compress(reactions, sortlist))

    return reactions


def check_reaction_input(reactant_list, reactive_substructs):
    """
    Checks for valid reaction input by checking the reactants
    contain the supplied substructures in the same order.

    Returns False if reaction inputs are incompatible with
    the reactive substructures provided.

    Parameters
    ----------
    reactant_list: NorthNet Compound Object

    reaction_template: NorthNet ReactionTemplate object

    Returns
    -------
    bool
        Whether the supplied list of compounds is compatible with the
        reaction substructure order supplied.
    """

    test_list = [
        r.Mol.HasSubstructMatch(s) for r, s in zip(reactant_list, reactive_substructs)
    ]

    return all(test_list)


def check_reaction_occurence(reactants, network, reaction_template):
    """
    Check if a reaction template has already been applied to a compound in a
    reaction network.

    Parameters
    ----------
    reactants: list or tuple of NorthNet Compound objects
    network: NorthNet Network object
    reaction_template: NorthNet ReactionTemplate object

    Returns
    -------
    bool
        Whether reaction type has been applied to the compound or not.
    """

    if not isinstance(reactants, tuple):
        # Written since Python was passing in the variable
        # inside tuples of length 1, rather than the tuple.
        reactants = (reactants,)

    reaction_name = reaction_template.Name

    result = []
    for reactant in reactants:
        used_reactions = reactant.Out
        used_classes = [network.get_reaction_name(r) for r in used_reactions]
        reaction_performed = reaction_name in used_classes
        result.append(reaction_performed)

    return all(result)
