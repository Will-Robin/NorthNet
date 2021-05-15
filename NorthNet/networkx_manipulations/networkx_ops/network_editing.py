def remove_nodes(G, node_list):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph from which nodes will be removed.
    node_list: list
        list of nodes to be removed from the graph.

    Returns
    -------
    F: networkx DiGraph
        Copy of G with nodes removed.

    '''
    F = G.copy()

    for n in node_list:
        if n in F.nodes:
            F.remove_node(n)

    return F

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
