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

def network_from_edge_list_csv(edges_file):
    '''
    Create a networkx DiGraph object using a list of edges from a .csv file.

    Parameters
    ----------
    edges_file: str
        Path to file containing edges
        (format: Source,Target,Type,Id,Label,timeset,Weight,4,5,6 newline)

    Returns
    -------
    G: networkx.DiGraph object
        network.
    '''

    # Build edges list
    edge_list = []
    with open(edges_file, "r") as f:
        line = next(f)
        for line in f:
            ln = line.strip("\n")
            ln = ln.split(",")
            edge_list.append(ln[:2])

    # Create networkx graph
    G = nx.DiGraph()

    # add edges to network.
    for e in edge_list:
        G.add_edge(e[0],e[1])

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    return G

def create_network_from_csv_files(coords_file, edges_file):
    '''
    Create a networkx DiGraph object using coordinates and edges files.

    Parameters
    ----------
    coords_file: str
        Path to file containing node coordinates (format: compound, x, y newline)
    edges_file: str
        Path to file containing edges
        (format: Source,Target,Type,Id,Label,timeset,Weight,4,5,6 newline)

    Returns
    -------
    G: networkx.DiGraph object
        network.
    '''

    G = network_view.network_from_edge_list_csv(edges_file)

    # Build coordinates list
    spec_coords = {}
    with open(coords_file, "r") as f:
        for line in f:
            ln = line.strip("\n")
            ln = ln.split(",")
            if ln[0] in G.nodes:
                pass
            else:
                G.add_node(ln[0])

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    G = network_view.set_network_coordinates(G, coords_file)

    return G
