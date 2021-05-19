def get_network_data(fname):
    print("Old function_name. Returning load_pickle_network_data()")
    return f_io.load_pickle_network_data(fname)

def load_pickle_network_data(fname):
    '''
    For loading pickled network data. Effectively a wrapper for pickle.
    Parameters
    ----------
    fname: str
        File to extract object from
    Returns
    -------
    net: object
        Object from pickle file.
    '''

    with open(fname, "rb") as f:
        net = pickle.load(f)

    return net

def read_gdf_file(file):
    '''
    Get network edges and coordinates from .gdf file
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

    with open(file, "r") as f:
        node_list = []

        for c,line in enumerate(f):
            if c == 0:
                spl = line.strip("\n").split(",")
                node_list_header = [x.split(" ")[-2] for x in spl]
            elif "edgedef" in line:
                break
            else:
                spl = line.strip("\n").split(",")
                node_list.append({k:v for k,v in zip(node_list_header, spl)})

    pos = {}
    for n in node_list:
        pos[n["name"]] = (float(n["x"]), float(n["y"]))

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

    edges = []
    for e in edge_list:
        edges.append((e[" node1"].strip("'"), e["node2"].strip("'")))

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

def read_gdf(file):

    with open(file, "r") as f:
        node_list = []

        for c,line in enumerate(f):
            if c == 0:
                spl = line.strip("\n").split(",")
                node_list_header = [x.split(" ")[-2] for x in spl]
            elif "edgedef" in line:
                break
            else:
                spl = line.strip("\n").split(",")
                node_list.append({k:v for k,v in zip(node_list_header, spl)})

    pos = {}
    for n in node_list:
        pos[n["name"]] = (float(n["x"]), float(n["y"]))

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

    edges = []
    for e in edge_list:
        edges.append((e[" node1"].strip("'"), e["node2"].strip("'")))

    return pos, edges

def get_species_coordinates(file):

    coords_dict = {}

    with open(file) as f:
        for line in f:
            ln = line.strip("\n")
            ln  = [x for x in ln.split(",") if x != ""]
            coords_dict[ln[0]+" M"] = [float(ln[1]),float(ln[2])]

    return coords_dict

def network_from_edge_list_csv(edges_file):
    '''
    Create a networkx DiGraph object using a list of edges from a .csv file.

    Parameters
    ----------
    edges_file: str
        Path to file containing edges
        (format: Source,Target,Type,Id,Label,timeset,Weight,4,5,6 newline)

    default_node_size: float

    default_node_color: float

    default_edge_width: int (float)
        default edge width
    default_edge_color: str
        Default edge color.

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

def coordinates_from_file(coords_file):
    '''
    Adds coordinates to node attributes from a .csv file.

    Parameters
    ----------
    G: networkx DiGraph
        Network
    coords_file: str
        Path to file containing node coordinates (format: compound, x, y newline)

    '''
    # Build coordinates list
    spec_coords = {}
    with open(coords_file, "r") as f:
        for line in f:
            ln = line.strip("\n")
            ln = ln.split(",")
            if ln[0].strip('"') in spec_coords:
                pass
            else:
                spec_coords[ln[0].strip('"')] = tuple([float(x) for x in ln[1:]])

    return spec_coords

def network_from_reaction_list(reaction_list, network_name = "", description = "", reaction_mapping = False):
    from NorthNet import Classes
    from rdkit import Chem
    from NorthNet import reaction_operations as r_ops
    from NorthNet import network_operations as net_ops

    header = ["Reaction", "Description", "Reaction ID","References"]
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    ID_dict = {S:str(c) for c,S in enumerate(alphabet)}
    network = Classes.Network([],network_name, description)
    add_to_net = []
    for r in reaction_list:

        reactants = [x for x in r.split(">>")[0].split(".") if x != ""]
        products = [x for x in r.split(">>")[1].split(".") if x != ""]

        r_mols = [Chem.MolFromSmiles(x) for x in reactants]
        p_mols = [Chem.MolFromSmiles(x) for x in products]

        r_canon = [Chem.MolToSmiles(x, isomericSmiles = True) for x in r_mols]
        p_canon = [Chem.MolToSmiles(x , isomericSmiles = True) for x in p_mols]

        r_string = r_ops.reconstitute_reaction(r_canon,p_canon)

        rID = "".join([ID_dict[s] for s in r_string if s in ID_dict])

        line = [r_string, "User input reaction.", rID, "None"]

        add_to_net.append(Classes.Reaction( Classes.Reaction_Database_Entry(header, line) ))

    network.add_reactions(add_to_net)

    if reaction_mapping == True:
        net_ops.map_network_reactions(network)

    return network
