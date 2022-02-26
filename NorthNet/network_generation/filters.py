def remove_reactions_by_product_substruct(network, substruct):
    """
    Remove the reactions in a nework if any of their products contain a defined
    substructure.

    Parameters
    ----------
    network: Classes.Network
        Network to check and modify.
    substruct: Classes.Substructure
        Substructure to check for in products.

    Returns
    -------
    None
    """

    remove_reactions = []

    for reaction in network.NetworkReactions:

        reaction_object = network.NetworkReactions[reaction]
        products = reaction_object.products

        for product in products:
            mol = network.NetworkCompounds[product].Mol

            if mol.HasSubstructMatch(substruct.Mol):
                remove_reactions.append(reaction_object)

    network.remove_reactions(remove_reactions)
