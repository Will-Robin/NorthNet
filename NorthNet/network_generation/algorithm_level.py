from rdkit import Chem
from NorthNet import Classes
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

    reactive_substructs = [Chem.MolFromSmarts(x)
                            for x in reaction_template.ReactantSubstructures]

    reactants = n_gen.reactive_species(list(network.NetworkCompounds.values()),
                                                            reactive_substructs)

    for r in reactants:
        insert = [r] + reagents
        if n_gen.check_reaction_input(insert, reactive_substructs):
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

    for r1 in reactants1:
        for r2 in reactants2:
            insert = [r1] + [r2]
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
    for x in reaction_template.ReactantSubstructures:
        substructure = Chem.MolFromSmarts(x)
        reactants.append(
            n_gen.reactive_species(
                        list(network.NetworkCompounds.values()), [substructure])
                        )

    # Build reactant combinations
    inputs = list(itertools.product(*reactants))
    for i in inputs:
        resulting_reactions = n_gen.run_reaction(i, reaction_template)
        network.add_reactions(resulting_reactions)
