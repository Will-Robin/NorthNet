def set_node_sizes(G, data_dict, scale_factor = 1e7, default_size = 1000.0, reaction_node_size = 250):
    '''
    To add size information into a networkx graph.

    CAUTION: will overwrite existing size information.

    Parameters
    ----------
    G: networkx DiGraph
        Network in which data will be inserted as node sizes.
    data_dict: dict
        Dictionary indexing node identities to a parameter to be used in node sizing.
    scale_factor: float
        Factor to mulitply data values by to make node size.
    default_size: int
        Size of nodes

    Returns
    -------
    G: networkx DiGraph
        Modified network.
    '''

    remove_units = {} # hack to convert keys from "X M" to "X"
    for d in data_dict:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_dict[d]

    for n in G.nodes:
        if n in remove_units:
            G.nodes[n]["size"] = remove_units[n]*scale_factor
        elif "S:" in n:
            G.nodes[n]["size"] = 800    #( G.degree(n, "width") - G.degree(n) )*20
        elif ">>" in n:
            G.nodes[n]["size"] = reaction_node_size
        else:
            G.nodes[n]["size"] = default_size

    return G

def set_node_borderwidths(G, data_dict, scale_factor = 10, default_size = 1.0):
    '''
    To add size information into a networkx graph.

    CAUTION: will overwrite existing size information.

    Parameters
    ----------
    G: networkx DiGraph
        Network in which data will be inserted as node sizes.
    data_dict: dict
        Dictionary indexing node identities to a parameter to be used in node sizing.
    scale_factor: float
        Factor to mulitply data values by to make node size.
    default_size: int
        Size of nodes

    Returns
    -------
    G: networkx DiGraph
        Modified network.
    '''

    remove_units = {} # hack to convert keys from "X M" to "X"
    for d in data_dict:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_dict[d]

    for n in G.nodes:
        if n in remove_units:
            G.nodes[n]["edgewidth"] = remove_units[n]*scale_factor
        elif "S:" in n:
            G.nodes[n]["edgewidth"] = default_size    #( G.degree(n, "width") - G.degree(n) )*20
        elif ">>" in n:
            G.nodes[n]["edgewidth"] = default_size
        else:
            G.nodes[n]["edgewidth"] = default_size

    return G


def set_node_colors(G, data_dict, default_color = "#3a393a"):
    '''
    Set the node colors according to node type.
    CAUTION: Overwites all previously added color information from network.

    Parameters
    ----------
    G: networkx DiGraph object
        Network on which to perform operation.

    data_dict: dict
        Dictionary with nodes as keys and data as values.

    default_color: str
        Color for nodes without assignments.

    Returns
    -------
    G: networkx DiGraph
        Modified graph with color attributes
    '''

    remove_units = {} # hack to convert keys from "X M" to "X"
    for d in data_dict:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_dict[d]

    for n in G.nodes:
        if n in remove_units:
            if n in info_params.colour_assignments:
                G.nodes[n]["color"] = info_params.colour_assignments[n]
            else:
                G.nodes[n]["color"] = "#000000"
        elif ">>" in n:
            G.nodes[n]["color"] = "#6a7a72"
        elif "S:" in n:
            G.nodes[n]["color"] = "r"
        else:
            G.nodes[n]["color"] = default_color

    return G

def set_node_shapes(G):

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["shape"]= "d"
        else:
            G.nodes[n]["shape"]= "o"

    return G

def set_node_colors_chain_length_gradient(G):
    '''
    Set the node colors according to node type.
    CAUTION: Overwites all previously added color information from network.

    Parameters
    ----------
    G: networkx DiGraph object
        Network on which to perform operation.

    Returns
    -------
    G: networkx DiGraph
        Modified graph with color attributes
    '''

    col_grad_dict = {0: 'k',
                     1: "#a7daf5",
                     2: '#75c9f5',
                     3: '#75a8f5',
                     4: '#4787e7',
                     5: '#425cea',
                     6: '#0520f8',
                     7: '#9b05f8'}

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["color"] = "#6a7a72"
        elif "S:" in n:
            G.nodes[n]["color"] = "r"
        else:
            G.nodes[n]["color"] = col_grad_dict[n.count("C")]

    return G

def set_edge_color(G, color = "#d1d1e0"):

    for e in G.edges:
        G.edges[e]["color"] = color

    return G

def set_edge_width(G, default_edge_width = 10):

    for e in G.edges:
        G.edges[e]["width"] = default_edge_width

    return G

def node_path_as_edge_widths(G, paths, rescale = 20):
    '''
    Parameters
    -----------
    G: networkx DiGraph
        Network
    paths: list
        list of paths to be superimpose on network
    rescale: float
        The maximum width of the edges.

    Returns
    -------
    Z: networkx DiGraph
        Network with paths superimposed.

    '''
    Z = G.copy() # copy the graph to make sure it is not modified

    for e in Z.edges: # assign a minimum weight attribute to all edges in network.
        Z.edges[e]["width"] = 1

    for p in paths:
        # iterate over paths. Every time an edge is traversed, 1 is added to the
        # width attribute
        for c,x in enumerate(p,1):
            if c == len(p):
                pass
            else:
                Z[x][p[c]]["width"] += 1
                Z[x][p[c]]["color"] = "r" # set a color so the edges pop out ot the eye.

    for e in Z.edges: # the edge widths are rescaled by the number of paths.
        Z.edges[e]["width"] /= len(paths)/rescale

    return Z

def edge_path_as_widths(G, paths, color = "#0925f0", rescale = 20):
    '''
    Parameters
    -----------
    G: networkx DiGraph
        Network
    paths: list
        list of paths to be superimpose on network
    rescale: float
        The maximum width of the edges.

    Returns
    -------
    Z: networkx DiGraph
        Network with paths superimposed.

    '''
    Z = G.copy() # copy the graph to make sure it is not modified

    for e in Z.edges: # assign a minimum weight attribute to all edges in network.
        Z.edges[e]["width"] = 5

    for p in paths:
        # iterate over paths. Every time an edge is traversed, 1 is added to the
        # width attribute
        for source,target in p:
            Z[source][target]["width"] += 1
            Z[source][target]["color"] = color # set a color so the edges pop out ot the eye.

    for e in Z.edges: # the edge widths are rescaled by the number of paths.
        Z.edges[e]["width"] /= len(paths)/rescale

    return Z

def apply_weight_distribution(graph, network, sigma, mu, attribute):
    '''
    To apply weights to graph edges according to reaction attributes.
    Graph is a Networkx graph, sigma is standard deviation, mu is the average,
    attribute is the numerical reaction attribute to be used for weighting.
    The weighting is according to a gaussian distribution with the mean set as
    the input reaction value. Since there are multiple reaction entries per
    reaction, the reactions are cycled through and the maximum weight calculated
    is selected for the edge weight.

    Parameters
    ----------
    graph: Networkx DiGraph object
        Graph to be weighted.

    network: NorthNet Network object
        NorthNet Network containing data for weighting.

    sigma: float
        Standard deviation of gaussian.
    mu: float
        Average of gaussian.
    attribute: str
        Attribute on which to base weighting.

    Returns
    -------
    graph: Networkx DiGraph object
        Graph with weights applied.
    '''

    for node in graph:
        if node in network.NetworkReactions:
            all_weights = []
            for x in network.NetworkReactions[node].Database_Entries:
                v = x.Info[attribute]
                all_weights.append(np.exp(-0.5*((v-mu)/sigma)**2))
            for r in network.NetworkReactions[node].Reactants:
                graph[r][node]["weight"] = max(all_weights)

    return graph
