def convert_to_networkx(network, default_node_size = 10, default_node_color = "b",
                        default_edge_width = 2, default_edge_color = "#6c6c7a", save_images = False):
    '''
    Converts NorthNet network object to networkx object.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be converted to networkx network.
    save_images: Bool
        Whether to create images or not.
    Returns
    -------
    G: networkx DiGraph object
        Networkx version of the NorthNet network.
    '''
    G = nx.DiGraph()
    n = 0
    C_patt = Chem.MolFromSmarts("[C]")
    for node in network.NetworkCompounds:
        if save_images:
            im = drw.create_image(node,n)
            n+=1
        else:
            im = "None"
        G.add_node( node, label = str(node), type = 'species', Image = im, C_count = len( network.NetworkCompounds[node].Mol.GetSubstructMatches(C_patt) ) )

    for r in network.NetworkReactions:

        ref_arr = []
        reaction_IDs = []
        for entry in network.NetworkReactions[r].Database_Entries:
            ref_arr.append(entry.Info["References"])
            reaction_IDs.append(entry.Info["Reaction ID"])

        for sr in network.NetworkReactions[r].Reactants:
            for sp in network.NetworkReactions[r].Products:
                G.add_edge(sr,r, Label = r, refs = ",".join(ref_arr), ReactionID = ",".join(reaction_IDs))
                G.add_edge(r,sp, Label = r, refs = ",".join(ref_arr), ReactionID = ",".join(reaction_IDs))

    for n in G.nodes:
        G.nodes[n]["size"] = default_node_size
        if ">>" in n:
            G.nodes[n]["color"] = "k"
        else:
            G.nodes[n]["color"] = default_node_color

    for e in G.edges:
        G.edges[e]["width"] = default_edge_width
        G.edges[e]["color"] = default_edge_color

    return G

def SNetwork_convert_to_networkx(snetwork, default_node_size = 20, default_node_color = "b",
                        default_edge_width = 2, default_edge_color = "#6c6c7a", save_images = False):
    '''
    Converts NorthNet network object to networkx object.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be converted to networkx network.
    save_images: Bool
        Whether to create images or not.
    Returns
    -------
    G: networkx DiGraph object
        Networkx version of the NorthNet network.
    '''
    G = nx.DiGraph()
    cont = 0
    C_patt = Chem.MolFromSmarts("[C]")
    for n in snetwork.SNetworkSubstructs:
        c_node_name  = "S:" + Chem.MolToInchi(snetwork.SNetworkSubstructs[n].Mol).replace(",",".")
        if save_images:
            im = drw.create_image(n, str(cont), fig_size = (150,150))
            cont+=1
        else:
            im = "None"

        G.add_node(c_node_name, label = n, type = 'substructure', Image = im)

    for n in snetwork.SNetworkCompounds:
        c_node_name  = snetwork.SNetworkCompounds[n].SMILES
        if save_images:
            im = drw.create_image(n, str(cont), fig_size = (150,150))
            cont+=1
        else:
            im = "None"

        G.add_node(c_node_name, label = n, type = "compound", Image = im, C_count = len( snetwork.SNetworkCompounds[n].Mol.GetSubstructMatches(C_patt) ))

        for i in snetwork.SNetworkCompounds[n].In:
            G.add_edge( "S:" + Chem.MolToInchi(snetwork.SNetworkSubstructs[i].Mol).replace(",","."), c_node_name)
        for o in snetwork.SNetworkCompounds[n].Out:
            G.add_edge(c_node_name, "S:" + Chem.MolToInchi(snetwork.SNetworkSubstructs[o].Mol).replace(",",".") )

    for r in snetwork.SNetworkTemplates:
        ref_arr = []
        reaction_IDs = []

        if len(snetwork.SNetworkTemplates[r].Database_Entries) == 0:
            r_node_name = snetwork.SNetworkTemplates[r].ReactionTemplate.Name
        else:
            r_node_name = snetwork.SNetworkTemplates[r].ReactionTemplate.Name

        G.add_node(r_node_name, label = r_node_name, type = "reaction")

        for entry in snetwork.SNetworkTemplates[r].Database_Entries:
            ref_arr.append(entry.Info["References"])
            reaction_IDs.append(entry.Info["Reaction ID"])

        for reac in snetwork.SNetworkTemplates[r].ReactionTemplate.ReactantSubstructures:
            G.add_edge("S:"+Chem.MolToInchi(snetwork.SNetworkSubstructs[reac].Mol).replace(",","."), r_node_name, Label = r_node_name, refs = ",".join(ref_arr), ReactionID = ",".join(reaction_IDs))

        for p in snetwork.SNetworkTemplates[r].ReactionTemplate.ProductSubstructures:
            G.add_edge(r_node_name, "S:"+Chem.MolToInchi(snetwork.SNetworkSubstructs[p].Mol).replace(",","."), Label = r_node_name, refs = ",".join(ref_arr), ReactionID = ",".join(reaction_IDs))

    for n in G.nodes:
        G.nodes[n]["size"] = default_node_size
        if ">>" in n:
            G.nodes[n]["color"] = "k"
        if "S:" in n:
            G.nodes[n]["color"] = "#da2034"
        else:
            G.nodes[n]["color"] = default_node_color

    for e in G.edges:
        G.edges[e]["width"] = default_edge_width
        G.edges[e]["color"] = default_edge_color

    return G


def simple_convert_to_networkx(network):
    print("simple_convert_to_networkx: changed name to convert_to_networkx. Returning that function instead.")
    return conv.convert_to_networkx(network, save_images = False)

def network_from_node_path(path, default_width = 5, default_color =  "#0925f0"):
    '''
    Create a networkx DiGraph object from a path (ordered list of nodes).

    Parameters
    ----------
    nodes: list
        list of nodes present in the network to be created (must have coordinates).

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
    # Build coordinates list
    G = nx.DiGraph()

    nx.add_path(G, path, width = default_width, color = default_color)

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    return G

def network_from_edge_path(path, default_width = 5, default_color =  "#0925f0"):
    '''
    Create a networkx DiGraph object from a path (list of edges).

    Parameters
    ----------
    nodes: list
        list of edges (tuples) present in the network to be created
        (must have coordinates).

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
    # Build coordinates list
    G = nx.DiGraph()

    G.add_edges_from(path, width = default_width, color = default_color)

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    return G

def network_from_multiple_node_paths(paths, default_width = 5, default_color =  "#0925f0"):
    '''
    Create a networkx DiGraph object from a path (ordered list of nodes).

    Parameters
    ----------
    paths: list of lists
        list containing lists of nodes present in the network to be created
        (must have coordinates).

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
    # Build coordinates list
    G = nx.DiGraph()

    for p in paths:
        G.add_edges_from(p, width = default_width, color = default_color)

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    return G

def network_from_multiple_edge_paths(paths, default_width = 5, default_color =  "#0925f0"):
    '''
    Create a networkx DiGraph object from a path (ordered list of nodes).

    Parameters
    ----------
    paths: list of lists
        list containing lists of edges (as tuples) present in the network to be
        created (must have coordinates).

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
    # Build coordinates list
    G = nx.DiGraph()

    for p in paths:
        nx.add_path(G, p, width = default_width, color = default_color)

    for n in G.nodes:
        if ">>" in n:
            G.nodes[n]["Type"] = "Reaction"
        elif "S:" in n:
            G.nodes[n]["Type"] = "Substructure"
        else:
            G.nodes[n]["Type"] = "Compound"

    return G

def node_path_to_edge_path(path, reaction_network):
    '''
    path: list
        A path in the form of a list of nodes.

    reaction_network: NorthNet Network object
        Reaction network holding reaction information.

    Returns
    --------
    edge_list: list
        List of edges (as tuples) including any secondary reactants required
        for the path (not including [OH-] or O).
    '''

    edge_list = []
    for c,p in enumerate(path,1):
        if c == len(path):
            pass
        elif ">>" in p:
            parent = reaction_network.NetworkReactions[p].ReactionSMILES

            reactant_edges = [(r, p) for r in reaction_network.NetworkReactions[p].Reactants if r != "[OH-]" and r != "O"]
            product_edges = [(p,prod) for prod in reaction_network.NetworkReactions[p].Products if prod != "[OH-]" and prod != "O"]
            edge_list.extend(reactant_edges)
            edge_list.extend(product_edges)
            edge_list.append((parent, path[c]))

        else:
            parent = p
            edge_list.append((parent, path[c]))

    return edge_list


def convert_path_to_substructure_path(path, substruct_network, reaction_network):
    '''

    Translates a "compound-reaction" path to a
    "compound-substructure-reaction_class-substructure-compound" path.

    Parameters
    ----------
    path: list (elements: str)
        A path of compounds and reactions.

    substruct_network: NorthNet Substructure_Network object
        Framework for translating the compound-reaction path to a
        compound-substructure-reaction_class-substructure-compound path.

    reaction_network: NorthNet Network object
        To help locate reaction templates.

    Returns
    -------
    edge_list: list (tuples)
        Path translation to edges in terms of substructure translations.
    '''

    edge_list = []
    for c,x in enumerate(path, 1):

        if c == len(path):
            # end of the sequence. Only need to append the reaction_template
            # the new path.
            pass

        elif ">>" in x:
            parent_name = reaction_network.NetworkReactions[x].ReactionTemplate.ReactionSMARTS
            possible_out_routes = reaction_network.NetworkReactions[x].ReactionTemplate.ProductSubstructures
            possible_linkers = substruct_network.SNetworkCompounds[path[c]].In
            linkers = list(set(possible_out_routes).intersection(possible_linkers))
            # substructures are labelled by their InChI names with commas
            # replaced with periods in networkx reprsentation.
            # This is a hack to deal with that issue.
            child_names = ["S:" + Chem.MolToInchi(substruct_network.SNetworkSubstructs[l].Mol).replace(",",".") for l in linkers]

            edge_set = [(parent_name, c_name) for c_name in child_names]

            edge_set.extend([(c_name,path[c]) for c_name in child_names])

            edge_list.extend(edge_set)

        else:
            parent_name = x
            # find which substucture leads to the next reaction
            possible_out_routes = substruct_network.SNetworkCompounds[x].Out
            possible_linkers = reaction_network.NetworkReactions[path[c]].ReactionTemplate.ReactantSubstructures

            linkers = list(set(possible_out_routes).intersection(possible_linkers))
            # substructures are labelled by their InChI names with commas
            # replaced with periods in networkx reprsentation.
            # This is a hack to deal with that issue.
            child_names = ["S:" + Chem.MolToInchi(substruct_network.SNetworkSubstructs[l].Mol).replace(",",".") for l in linkers]
            edge_set = [(parent_name, c_name) for c_name in child_names]

            edge_set.extend([(c_name,reaction_network.NetworkReactions[path[c]].ReactionTemplate.ReactionSMARTS) for c_name in child_names])

            edge_list.extend(edge_set)

    return edge_list


def convert_path_to_reaction_class_path(path, substruct_network, reaction_network):
    '''

    Translates a "compound-reaction" path to a
    "compound-substructure-reaction_class-substructure-compound" path.

    Parameters
    ----------
    path: list (elements: str)
        A path of compounds and reactions.

    substruct_network: NorthNet Substructure_Network object
        Framework for translating the compound-reaction path to a
        compound-substructure-reaction_class-substructure-compound path.

    reaction_network: NorthNet Network object
        To help locate reaction templates.

    Returns
    -------
    edge_list: list (tuples)
        Path translation to edges in terms of substructure translations.
    '''

    edge_list = []

    for c,x in enumerate(path, 1):

        if c == len(path):
            # end of the sequence. Only need to append the reaction_template
            # the new path.
            pass

        elif ">>" in x:
            parent_name = reaction_network.NetworkReactions[x].ReactionTemplate.ReactionSMARTS
            out_routes = reaction_network.NetworkReactions[x].Products
            edges = [(parent_name, o) for o in out_routes]
            [edge_list.append(e) for e in edges]

        else:
            parent_name = x
            # find which substucture leads to the next reaction
            next_reaction = reaction_network.NetworkReactions[path[c]].ReactionTemplate.ReactionSMARTS
            edge = [(parent_name, next_reaction)]
            [edge_list.append(e) for e in edge]


    return edge_list
