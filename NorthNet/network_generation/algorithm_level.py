from rdkit import Chem
from NorthNet import network_generation as n_gen

def extend_network_specific(network, reagents, reaction_template):
    '''
    Extend the network using a single reagent set.

    Parameters
    ----------
    network: NorthNet Network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.

    Returns
    -------
    None
    '''

    template_substructures = reaction_template.ReactantSubstructures
    compounds_in_network = list(network.NetworkCompounds.values())
     
    reactive_substrs = [Chem.MolFromSmarts(x) for x in template_substructures]

    reactants = n_gen.reactive_species(compounds_in_network, reactive_substrs)

    for reactant in reactants:
        insert = [reactant] + reagents
        reaction_done = n_gen.check_reaction_occurence(
                                                        reactant, 
                                                        network, 
                                                        reaction_template
                                                        )
        input_valid = n_gen.check_reaction_input(insert, reactive_substrs)

        if not reaction_done and input_valid:
            resulting_reactions = n_gen.run_reaction(insert, reaction_template)
            network.add_reactions(resulting_reactions)
        else:
            pass

def extend_network_self(network, reaction_template):
    '''
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: NorthNet Network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    secondary_substructure: NorthNet/NetGen Substructure object
        Substructure of the second reaction component.

    Returns
    -------
    None
    '''
    substructures = [Chem.MolFromSmarts(x)
                            for x in reaction_template.ReactantSubstructures]

    reactants1 = n_gen.reactive_species(
                    list(network.NetworkCompounds.values()), [substructures[0]]
                    )

    reactants2 = n_gen.reactive_species(
                    list(network.NetworkCompounds.values()), [substructures[1]]
                    )

    for reactant_1 in reactants1:
        reaction_done = n_gen.check_reaction_occurence(
                                                        reactant_1, 
                                                        network, 
                                                        reaction_template
                                                        )
        if reaction_done:
            continue

        for reactant_2 in reactants2:
            insert = [reactant_1] + [reactant_2]
            resulting_reactions = n_gen.run_reaction(insert, reaction_template)
            network.add_reactions(resulting_reactions)

def extend_network_task(network, reaction_template):
    '''
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: NorthNet Network object
        Network to be extrapolated from. Modified in place.
    reagents: NorthNet Compound objects
        Reagents to be applied to the network.
    reaction_template: NorthNet Reaction_Template object.
        Reaction template to be used on the network.
    secondary_substructure: NorthNet Substructure object
        Substructure of the second reaction component.

    Returns
    -------
    None
    '''
    import itertools

    reactants = []
    for substruct in reaction_template.ReactantSubstructures:
        substructure = Chem.MolFromSmarts(substruct)
        reactants.append(
            n_gen.reactive_species(
                        list(network.NetworkCompounds.values()), [substructure])
                        )

    # Build reactant combinations
    inputs = list(itertools.product(*reactants))
    for input in inputs:
        reaction_done = n_gen.check_reaction_occurence(
                                                        input[0], 
                                                        network, 
                                                        reaction_template
                                                        )
        if not reaction_done:
            resulting_reactions = n_gen.run_reaction(input, reaction_template)
            network.add_reactions(resulting_reactions)

