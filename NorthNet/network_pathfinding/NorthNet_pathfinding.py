from NorthNet import network_pathfinding as pathf
import networkx as nx
from NorthNet import Classes
from NorthNet import conversions as conv

def next_node_layer(node,network):
    '''
    Takes a node in a network object, blocks it and finds the products of
    reactions leading from it using the 'Out' edges. It then creates a
    list of unblocked nodes, transferring the reactions logged in the original
    node into them along with the new reactions.
    The reaction log of the input node is cleared and a list
    of new nodes is returned.
    Parameters
    ----------
    node: str
        SMILES of node to start expansion from.

    network: NorthNet Network object
        NorthNet Network to be searched.

    Returns
    -------
    nextnodes: list
        List of SMILES structures of the next node layer.
    '''

    network.NetworkCompounds[node].Block = True
    nextnodes = []
    edges = network.NetworkCompounds[node].Out

    for e in edges:
        prods = network.NetworkReactions[e].Products
        if len(prods) < 100 and len(network.NetworkReactions[e].Reactants) < 100:
            for pr in prods:
                if network.NetworkCompounds[pr].Block:
                    pass
                elif abs(pr.count("C") + pr.count("c")-node.count("C")-node.count("c")) > 1:
                    pass
                else:
                    '''Add the reactions to a node EdgeLog'''
                    network.NetworkCompounds[pr].EdgeLog += network.NetworkCompounds[node].EdgeLog
                    network.NetworkCompounds[pr].EdgeLog.append(e)
                    network.NetworkCompounds[pr].EdgeLog = list(set(network.NetworkCompounds[pr].EdgeLog))
                    nextnodes.append(pr)

    nextnodes = list(set(nextnodes)) # remove duplicates
    #network.NetworkCompounds[node]['EdgeLog'] = [] # clear previous edgelog

    return nextnodes

def next_node_layer_2(root_node, network):
    '''
    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be searched.
    Returns
    -------
    nextnodes: list
        List of SMILES structures of the next node layer.
    '''
    nextnodes = set()

    if network.NetworkCompounds[root_node].Block == True:
        pass
    else:
        network.NetworkCompounds[root_node].Block = True

        outgoing = network.NetworkCompounds[root_node].Out

        for o in outgoing:
            for p in network.NetworkReactions[o].Products:
                if network.NetworkCompounds[p].Block == True:
                    pass
                else:
                    nextnodes.add(p)
                    '''Add the reactions to a node EdgeLog'''
                    network.NetworkCompounds[p].EdgeLog += network.NetworkCompounds[root_node].EdgeLog
                    network.NetworkCompounds[p].EdgeLog.append(entry)
                    network.NetworkCompounds[p].EdgeLog = list(set(network.NetworkCompounds[p].EdgeLog))

    network.NetworkCompounds[root_node].EdgeLog = []

    return nextnodes

def expand_around_node(network,root_nodes, expansion_steps = 1):
    '''Expands around a root_node in a network object for n expansion_steps using
    the next_node_layer function. A current_layer is created and overwritten
    for each iteration. If the root_node is found in the current_layer, a cycle
    has been found and the function returns the reactions logged in the
    root_node. If a cycle is not found the function continues until all of the
    iterations have been completed (or memory runs out). In this case the
    function returns an empty list.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be searched.
    root_nodes: list
        List of reaction SMILES for nodes to start from.
    expansion_steps: int
        Number of expansion iterations to perform.
    Returns
    -------
    new_edges: list:
        List of reactions.
    '''

    current_layer = root_nodes
    new_edges = []
    n = 0
    while n < expansion_steps:
        next_layer = []

        for l in current_layer:
            nl = pathf.next_node_layer(l, network)
            next_layer.extend(nl)

        current_layer = next_layer

        if n == 4:
            for node in root_nodes:
                network.NetworkCompounds[node].Block = False

        if any(item in current_layer for item in root_nodes):

            print('found cycle?')
            sub_key_list = network.NetworkCompounds[root_nodes[0]].EdgeLog
            sub_key_list = list(set(sub_key_list))
            print('Number of edges: ',len(sub_key_list))

            for node in network.NetworkCompounds:
                network.NetworkCompounds[node].Block = False

            new_edges = []
            [new_edges.extend(network.NetworkReactions[k]) for k in sub_key_list]
            break

        n += 1

    return new_edges

def expand_around_node_2(network,root_nodes, expansion_steps = 1):
    '''
    Expands around a root_node in a network object for n expansion_steps using
    the next_node_layer function. A current_layer is created and overwritten
    for each iteration. If the root_node is found in the current_layer, a cycle
    has been found and the function returns the reactions logged in the
    root_node. If a cycle is not found the function continues until all of the
    iterations have been completed (or memory runs out). In this case the
    function returns an empty list.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be searched.
    root_nodes: list
        List of reaction SMILES for nodes to start from.
    expansion_steps: int
        Number of expansion iterations to perform.
    Returns
    -------
    new_edges: list:
        List of reactions.
    '''

    current_layer = root_nodes
    new_edges = []
    n = 0
    while n < expansion_steps:

        print(n)
        next_layer = []

        for l in current_layer:
            nl = pathf.next_node_layer_2(l, network)
            next_layer.extend(nl)

        current_layer = next_layer

        if n == 2:
            for node in root_nodes:
                network.NetworkCompounds[node].Block = False

        if any(item in current_layer for item in root_nodes):

            print('found cycle?')
            reaction_cycle = network.NetworkCompounds[root_nodes[0]].EdgeLog
            print('Number of edges: ',len(reaction_cycle))

            return Classes.Network(reaction_cycle)
        n += 1

    return False

def pair_cycle_search(network,root_node):
    '''
    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be searched.
    root_node: str
        SMILES for node to start from.
    Returns
    -------
    None
    '''
    G = conv.create_networkx_graph(network)

    sec_layer = []
    for o in network.NetworkCompounds[root_node].Out:
        for p in network.NetworkReactions[o].Products:
            sec_layer.append(p)

    for s in sec_layer:
        for pt in nx.all_simple_paths(G, source = root_node, target = s):
            print(pt)

def new_cycle_finder(network, root_nodes, expansion_steps = 10):

    '''
    Expands around a root_node in a network object for n expansion_steps using
    the next_node_layer function. A current_layer is created and overwritten
    for each iteration. If the root_node is found in the current_layer, a cycle
    has been found and the function returns the reactions logged in the
    root_node. If a cycle is not found the function continues until all of the
    iterations have been completed (or memory runs out). In this case the
    function returns an empty list.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be searched.
    root_node: str
        SMILES for node to start from.
    expansion_steps: int
        Number of expansion iterations to perform.
    Returns
    -------
    output: list
        List of reactions.
    '''

    current_layer = root_nodes
    new_edges = []
    n = 0
    while n < expansion_steps:
        print(n)
        next_layer = []
        reactions = []

        for l in current_layer:
            if network.NetworkCompounds[l].Block:
                reactions = network.NetworkCompounds[l].EdgeLog
                n = expansion_steps
                continue
            else:
                network.NetworkCompounds[node].Block = True
                edges = network.NetworkCompounds[node].Out

                for e in edges:
                    prods = network.NetworkReactions[e].Products
                    for pr in prods:
                        '''Add the reactions to a node EdgeLog'''
                        nextnodes.append(pr)

                nextnodes = list(set(nextnodes)) # remove duplicates
                #network.NetworkCompounds[node]['EdgeLog'] = [] # clear previous edgelog
                next_layer.extend(nl)

        current_layer = next_layer
        n += 1

    output = [network.NetworkReactions[x] for x in reactions]

    return output
