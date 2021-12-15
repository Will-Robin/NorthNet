import json
from graphviz import Digraph

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

    # Create a graph with graphviz to plot a scheme of the network
    dot = Digraph(comment = '',
                  engine = render_engine,
                  strict = 'True',
                  format = 'json')

    for node in network.NetworkCompounds:
        dot.node(node,node)

    for reaction in network.NetworkReactions: # Create edges between the nodes.

        reactants = network.NetworkReactions[reaction].Reactants
        products = network.NetworkReactions[reaction].Products

        dot.node(reaction, reaction)

        for reactant in reactants:
            dot.edge(reactant,reaction)

        for product in products:
            dot.edge(reaction,product)

    json_string = dot.pipe().decode()

    layout_in_string = json.loads(json_string)

    pos = {}
    for l_obj in layout_in_string['objects']:
        pos[l_obj['name']] = [float(x) for x in l_obj['pos'].split(',')]

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

    # Create a graph with graphviz to plot a scheme of the network
    dot = Digraph(comment = '', engine=render_engine, strict = 'True',
                format = 'json')

    for node in network.nodes: # add in nodes
        dot.node(node,node)

    for edge in network.edges: # Create edges between the nodes.
        dot.edge(edge[0],edge[1])

    json_string = dot.pipe().decode()

    layout_in_string = json.loads(json_string)

    pos = {}
    for l_obj in layout_in_string['objects']:
        pos[l_obj['name']] = [float(x) for x in l_obj['pos'].split(',')]

    return pos
