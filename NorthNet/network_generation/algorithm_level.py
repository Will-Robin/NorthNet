from rdkit import Chem
from NorthNet import network_generation as n_gen
from NorthNet import Classes

def extend_network_specific(network, reagents, reaction_template, exceptions):
    '''
    Extend the network using a single reagent set.

    Parameters
    ----------
    network: NorthNet network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    exceptions: list
        List of forbidden substructures.

    Returns
    -------
    None
    '''

    reactive_substructs = [Chem.MolFromSmarts(x) for x in reaction_template.ReactantSubstructures]
    reactants = n_gen.reactive_species(list(network.NetworkCompounds.values()), reactive_substructs)

    for r in reactants:
        insert = [r] + reagents
        if n_gen.check_reaction_input(insert, reactive_substructs):
            deprot = n_gen.run_reaction(insert, reaction_template)
            deprot = n_gen.remove_invalid_reactions(deprot, exceptions)
            insertion = [Classes.Reaction(r) for r in deprot]
            network.add_reactions(insertion) # extend the reactions list with the new reactions
        else:
            pass

def intra_network_reactions(network,  reagents, reaction_template, exceptions):
    '''
    Extend the network using a single reagent set.

    Parameters
    ----------
    network: NorthNet network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    exceptions: list
        List of forbidden substructures.

    Returns
    -------
    None
    '''

    reactive_substructs = [Chem.MolFromSmarts(x) for x in reaction_template.ReactantSubstructures]
    reactants = n_gen.reactive_species(list(network.NetworkCompounds.values()), reactive_substructs)

    reactions_to_add = []
    for r in reactants:
        insert = [r] + reagents
        if n_gen.check_reaction_input(insert, reactive_substructs):
            deprot = n_gen.run_reaction(insert, reaction_template)
            deprot = n_gen.remove_invalid_reactions(deprot, exceptions)
            reactions_to_add.extend([Classes.Reaction(r) for r in deprot])
        else:
            pass

    return reactions_to_add

def extend_network_self(network, reaction_template, exceptions):
    '''
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: NorthNet network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    secondary_substructure: NorthNet/NetGen Substructure object
        Substructure of the second reaction component.
    exceptions: list
        List of forbidden substructures.

    Returns
    -------
    None
    '''
    substructures = [Chem.MolFromSmarts(x) for x in reaction_template.ReactantSubstructures]

    reactants1 = n_gen.reactive_species(list(network.NetworkCompounds.values()), [substructures[0]]) # cruder way of getting the SMARTS back
    reactants2 = n_gen.reactive_species(list(network.NetworkCompounds.values()), [substructures[1]])

    for r1 in reactants1:
        for r2 in reactants2:
            insert = [r1] + [r2]
            deprot = n_gen.run_reaction(insert, reaction_template)
            deprot = n_gen.remove_invalid_reactions(deprot, exceptions)
            insertion = [Classes.Reaction(r) for r in deprot]
            network.add_reactions(insertion)

def extend_network_self_with_constraint(network, reaction_template, exceptions):
    '''
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: NorthNet network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    secondary_substructure: NetGen Substructure object
        Substructure of the second reaction component.
    exceptions: list
        List of forbidden substructures.

    Returns
    -------
    None
    '''
    substructures = [Chem.MolFromSmarts(x) for x in reaction_template.ReactantSubstructures]

    reactants1 = n_gen.reactive_species(list(network.NetworkCompounds.values()), [substructures[0]]) # cruder way of getting the SMARTS back
    reactants2 = n_gen.reactive_species(list(network.NetworkCompounds.values()), [substructures[1]])

    C_patt = Chem.MolFromSmarts('[C]')
    reactants1 = [r for r in reactants1 if len(r.Mol.GetSubstructMatch(C_patt)) < 4]
    reactants2 = [r for r in reactants2 if len(r.Mol.GetSubstructMatch(C_patt)) < 4]

    for r1 in reactants1:
        for r2 in reactants2:
            insert = [r1] + [r2]
            deprot = n_gen.run_reaction(insert, reaction_template)
            deprot = n_gen.remove_invalid_reactions(deprot, exceptions)
            insertion = [Classes.Reaction(r) for r in deprot]
            network.add_reactions(insertion)

def extend_network_task(network, reaction_template, exceptions):
    '''
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: NorthNet network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    secondary_substructure: NorthNet Substructure object
        Substructure of the second reaction component.
    exceptions: list
        List of forbidden substructures.

    Returns
    -------
    None
    '''
    import itertools

    substructures = [Chem.MolFromSmarts(x) for x in reaction_template.ReactantSubstructures]

    reactant_num = len(substructures)

    reactants = []
    for s in substructures:
        reactants.append(n_gen.reactive_species(list(network.NetworkCompounds.values()), [s]))

    # Build reactant combinations
    inputs = list(itertools.product(*reactants))
    for i in inputs:
        deprot = n_gen.run_reaction(i, reaction_template)
        deprot = n_gen.remove_invalid_reactions(deprot, exceptions)
        insertion = [Classes.Reaction(r) for r in deprot]
        network.add_reactions(insertion)
