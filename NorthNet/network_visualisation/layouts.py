def graphviz_layout_NorthNet(network, render_engine = 'fdp'):
    '''
    Uses graphviz to generate a layout from a NorthNet Network.

    Layouts from the graphviz documentation:
    dot − filter for drawing directed graphs
    neato − filter for drawing undirected graphs
    twopi − filter for radial layouts of graphs
    circo − filter for circular layout of graphs
    fdp − filter for drawing undirected graphs
    sfdp − filter for drawing large undirected graphs
    patchwork − filter for squarified tree maps
    osage − filter for array-based layouts

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be represented.
    render_engine: str
        Layout render engine for graphviz.
    Returns
    -------
    pos: dict
        Compounds SMILES and Reaction SMILES are keys to their positions.
        {SMILES:[float(x),float(y)]}
    '''
    import json
    from graphviz import Digraph

    # Create a graph with graphviz to plot a scheme of the network
    dot = Digraph(comment = '',
                  engine = render_engine,
                  strict = 'True',
                  format = 'json')

    for n in network.NetworkCompounds:
        dot.node(n,n)

    for r in network.NetworkReactions: # Create edges between the nodes.

        reactants = network.NetworkReactions[r].Reactants
        products = network.NetworkReactions[r].Products

        dot.node(r, r)

        for reac in network.NetworkReactions[r].Reactants:
            dot.edge(reac,r)

        for p in network.NetworkReactions[r].Products:
            dot.edge(r,p)

    json_string = dot.pipe().decode()

    y = json.loads(json_string)

    pos = {}
    for o in y['objects']:
        pos[o['name']] = [float(x) for x in o['pos'].split(',')]

    return pos

def graphviz_layout_networkx(network, render_engine = 'fdp'):
    '''
    Uses graphviz to generate a layout from a networkx graph.

    Layouts from the graphviz documentation:
    dot − filter for drawing directed graphs
    neato − filter for drawing undirected graphs
    twopi − filter for radial layouts of graphs
    circo − filter for circular layout of graphs
    fdp − filter for drawing undirected graphs
    sfdp − filter for drawing large undirected graphs
    patchwork − filter for squarified tree maps
    osage − filter for array-based layouts

    Parameters
    ----------
    network: Networkx DiGraph
        Network to be represented.
    render_engine: str
        Layout render engine for graphviz.
    Returns
    -------
    pos: dict
        Compounds SMILES and Reaction SMILES are keys to their positions.
        {SMILES:[float(x),float(y)]}
    '''
    import json
    from graphviz import Digraph

    # Create a graph with graphviz to plot a scheme of the network
    dot = Digraph(comment = '', engine=render_engine, strict = 'True',
                format = 'json')

    for n in network.nodes: # add in nodes
        dot.node(n,n)

    for e in network.edges: # Create edges between the nodes.
        dot.edge(e[0],e[1])

    json_string = dot.pipe().decode()

    y = json.loads(json_string)

    pos = {}
    for o in y['objects']:
        pos[o['name']] = [float(x) for x in o['pos'].split(',')]

    return pos
