def load_network_from_reaction_list(reaction_list, name = '', description = ''):
    '''
    Create a NorthNet Network from a list of reactions

    reaction_list: list of reaction SMILES strings
        Format: e.g. C=O.OC=C(O)CO>>O=C([C@@H](CO)O)CO

    network: NorthNet Network object
    '''
    from NorthNet import Classes

    rxns = []
    for r in reaction_list:
        rxns.append(Classes.Reaction(r))

    network = Classes.Network(rxns, name, description)

    return network

def node_edge_list_from_gdf(file):
    '''
    Get network edges and coordinates from .gdf file

    file: str or pathlib Path
        Path to file

    node_list, edge_list: list of dicts

    '''
    with open(file, "r") as f:
        coordinates = []

        for c,line in enumerate(f):
            if c == 0:
                spl = line.strip("\n").split(",")
                node_list_header = [x.split(" ")[-2] for x in spl]
            elif "edgedef" in line:
                break
            else:
                spl = line.strip("\n").split(",")
                coordinates.append({k:v for k,v in zip(node_list_header, spl)})

    node_list = {k["name"].strip(" "):[float(k["x"]),float(k["y"])] for k in coordinates}

    with open(file, "r") as f:
        edge_list = []
        readstate = False
        for c,line in enumerate(f):
            if "edgedef" in line:
                spl = line.strip("\n").strip(">edgedef").split(",")
                edge_list_header = [x for x in spl]

                readstate = True
            elif readstate:
                spl = line.strip("\n").split(",")
                edge_list.append({k:v for k,v  in zip(edge_list_header, spl)})

    return node_list, edge_list

def network_from_gdf(file):

    nodes, edges = node_edge_list_from_gdf(file)

    # Create networkx graph
    G = nx.DiGraph()

    # add edges to network.
    for e in edges:
        G.add_edge(e[0],e[1])

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    return G
