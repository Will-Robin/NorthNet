def remove_long_reactions(graph,network, length = 10):
    '''
    Parameters
    ----------
    graph: Networkx DiGraph object
        Graph to be mnodified.

    network: NorthNet Network object
        NorthNet Network containing reactions.
    length: int
        Maxium length of reactions to retain.
    Returns
    -------
    graph: Networkx DiGraph object
        Graph with long reactions removed.
    '''
    remove_list = []
    for r in network.NetworkReactions:
        if len(network.NetworkReactions[r].Products) > length:
            remove_list.append(r)
        elif len(network.NetworkReactions[r].Reactants) > length:
            remove_list.append(r)

    for rem in remove_list:
        graph.remove_node(rem)

    return graph
