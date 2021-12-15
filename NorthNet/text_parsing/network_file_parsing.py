from NorthNet import Classes

def load_network_from_reaction_list(reaction_list, name = '', description = ''):
    '''
    Create a NorthNet Network from a list of reactions

    reaction_list: list of reaction SMILES strings
        Format: e.g. C=O.OC=C(O)CO>>O=C([C@@H](CO)O)CO

    network: NorthNet Network object
    '''

    rxns = []
    for reaction in reaction_list:
        rxns.append(Classes.Reaction(reaction))

    network = Classes.Network(rxns, name, description)

    return network

def node_edge_list_from_gdf(filename):
    '''
    Get network edges and coordinates from .gdf file

    file: str or pathlib Path
        Path to file

    node_list, edge_list: list of dicts

    '''
    with open(filename, "r") as file:
        coordinates = []
        node_list_header = []

        for c,line in enumerate(file):
            if c == 0:
                spl = line.strip("\n").split(",")
                node_list_header = [x.split(" ")[-2] for x in spl]
                continue

            if "edgedef" in line:
                break

            spl = line.strip("\n").split(",")
            coordinates.append({k:v for k,v in zip(node_list_header, spl)})

    node_list = {k["name"].strip(" "):[float(k["x"]),float(k["y"])] for k in coordinates}

    with open(filename, "r") as file:
        edge_list = []
        edge_list_header = []

        readstate = False
        for c,line in enumerate(file):
            if "edgedef" in line:
                spl = line.strip("\n").strip(">edgedef").split(",")
                edge_list_header = spl

            if readstate:
                spl = line.strip("\n").split(",")
                edge_list.append({k:v for k,v  in zip(edge_list_header, spl)})

    return node_list, edge_list

def network_from_gdf(filename):

    import networkx as nx

    _, edges = node_edge_list_from_gdf(filename)


    # Create networkx graph
    G = nx.DiGraph()

    # add edges to network.
    for edge in edges:
        G.add_edge(edge [0],edge[1])

    for node in G.nodes:
        if ">>" in node:
            G.nodes[node]["Type"] = "Reaction"
        elif "S:" in node:
            G.nodes[node]["Type"] = "Substructure"
        else:
            G.nodes[node]["Type"] = "Compound"

    return G
