def paths_from_parameters_just_from_root(G, data_container, root_node  = "O=C(CO)CO"):
    '''
    Parameters
    ----------
    G: NetworkX DiGraph
        Graph in which paths will be found.
    data_container: dict
        Dictionary containing npode data.

    Returns
    -------
    paths: list of paths
    '''

    remove_units = {} # to convert keys frpom "X M" to "X"
    for d in data_container:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_container[d]

    # order the data in the dictionary by values
    ordered_data = [x for v,x in sorted( zip( list(remove_units.values()), [*remove_units] ) ) if x != "C=O"]
    ordered_data.reverse() # ordered in decreasing parameter value.
    ordered_data.remove(root_node) # remove the root node so it can be reinserted at the top of the list as a source.
    ordered_data.insert(0,root_node)

    Z = G.copy() # create a copy of the network
    Z.remove_node("C=O") # remove formaldehyde so that not all pathways pass through it.

    remove_list = []
    for e in Z.edges:
        for x in e:
            if ">>" in x:
                if "C=O" in x.split(">>")[1]:
                    remove_list.append(e)

    for e in remove_list:
        Z.remove_edge(e[0],e[1])

    paths = []
    for x in range(0,len(ordered_data)):
        if ordered_data[x] not in Z.nodes:
            pass
        elif nx.has_path(Z, root_node, ordered_data[x]): # test if there is a path between the compounds
            # Link together consecutive compounds in ordered_data
            path = nx.shortest_path(Z, source = root_node, target = ordered_data[x])
            paths.append(path)
        else: # if there is no path, pass
            pass

    return paths

def paths_from_parameters(G, data_container, root_node  = "O=C(CO)CO"):
    '''
    General idea: find shortest path to a node from any node upstream from it
    (when nodes are ordered in decreasing amplitude).

    Parameters
    ----------
    G: NetworkX DiGraph
        Graph in which paths will be found.
    data_container: dict
        Dictionary containing node data.
    root_node: str
        The compound from which the network grows.

    Returns
    -------
    paths: list of paths
    '''

    remove_units = {} # to convert keys frpom "X M" to "X"
    for d in data_container:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_container[d]

    # order the data in the dictionary by values (low to high)
    ordered_data = [x for v,x in sorted( zip( list(remove_units.values()), [*remove_units] ) ) if x != "C=O"]
    ordered_data.reverse() # reverse so order is high to low.
    ordered_data.remove(root_node) # remove the root node so it can be rsinserted at the top of the list as a source.
    ordered_data.insert(0,root_node) # place the root node at the top of the list

    Z = G.copy() # create a copy of the network
    Z.remove_node("C=O") # remove formaldehyde so that not all pathways pass through it.

    path_dict = {} # to store paths
    for x in range( 0, len(ordered_data) - 1 ):
        # x will be the source node
        for y in range( x+1, len(ordered_data)):
            # y is the target node and is selected from all the nodes downstream from x
            if ordered_data[y] not in path_dict:
                # create an entry in which to store paths
                path_dict[ordered_data[y]] = []

            if nx.has_path(Z, ordered_data[x], ordered_data[y]): # test if there is a path between the compounds
                # Link together consecutive compounds in ordered_data with shortest path
                path = nx.shortest_path(Z, source = ordered_data[x], target = ordered_data[y])
                path_dict[ordered_data[y]].append(path) # store

    paths = [] # to store output paths

    for p in path_dict: # finding the shortest path leading to the node. If there is more than one, one is taken and the others are ignored.
        paths.append(min(path_dict[p], key = len))

    return paths

def paths_from_parameters_min_amps_score(G, data_container, root_node  = "O=C(CO)CO"):
    '''
    General idea: find shortest path to a node from any node upstream from it
    (when nodes are ordered in decreasing amplitude).

    Parameters
    ----------
    G: NetworkX DiGraph
        Graph in which paths will be found.
    data_container: dict
        Dictionary containing node data.
    root_node: str
        The compound from which the network grows.

    Returns
    -------
    paths: list of paths
    '''

    remove_units = {} # to convert keys frpom "X M" to "X"
    for d in data_container:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_container[d]

    # order the data in the dictionary by values (low to high)
    ordered_data = [x for v,x in sorted( zip( list(remove_units.values()), [*remove_units] ) ) if x != "C=O"]
    ordered_data.reverse() # reverse so order is high to low.
    ordered_data.remove(root_node) # remove the root node so it can be rsinserted at the top of the list as a source.
    ordered_data.insert(0,root_node) # place the root node at the top of the list

    Z = G.copy() # create a copy of the network
    Z.remove_node("C=O") # remove formaldehyde so that not all pathways pass through it.

    path_dict = {} # to store paths
    for x in range( 0, len(ordered_data) - 1 ):
        # x will be the source node
        for y in range( x+1, len(ordered_data)):
            # y is the target node and is selected from all the nodes downstream from x
            if ordered_data[y] not in path_dict:
                # create an entry in which to store paths
                path_dict[ordered_data[y]] = []

            if nx.has_path(Z, ordered_data[x], ordered_data[y]): # test if there is a path between the compounds
                # Link together consecutive compounds in ordered_data with shortest path
                path = nx.shortest_path(Z, source = ordered_data[x], target = ordered_data[y])
                path_dict[ordered_data[y]].append(path) # store

    paths = [] # to store output paths

    for p in path_dict: # finding the shortest path leading to the node. If there is more than one, one is taken and the others are ignored.
        scores = []
        for poss in path_dict[p]:
            score = 0
            for x in poss:
                if ">>" in x:
                    pass
                elif x in remove_units:
                    score += remove_units[x]
                else:
                    score += 0
            scores.append(score)

        ind = scores.index(max(scores))
        paths.append(path_dict[p][ind])

    return paths

def paths_from_parameters_decay_score(G, data_container, motifs = False, root_node  = "O=C(CO)CO", depth = 20):
    '''
    General idea: find shortest path to a node from any node upstream from it
    (when nodes are ordered in decreasing amplitude).

    Parameters
    ----------
    G: NetworkX DiGraph
        Graph in which paths will be found.
    data_container: dict
        Dictionary containing node data.
    motifs: dict or False
        {"name":{"connection_point":[], "graph": nx.DiGraph}}
        Motifs are removed for the pathway search process
    root_node: str
        The compound from which the network grows.

    Returns
    -------
    paths: list of paths
    '''

    remove_units = {} # to convert keys frpom "X M" to "X"
    for d in data_container:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_container[d]

    remove_units[root_node] = 10.0 # force root_node to be highest amplitude

    Z = G.copy() # create a copy of the network
    if 'C=O' in Z.nodes:
        Z.remove_node("C=O") # remove formaldehyde so that not all pathways pass through it.

    # order the data in the dictionary by values (low to high)
    ordered_data = [x for v,x in sorted( zip( list(remove_units.values()), [*remove_units] ) ) if x != "C=O"]
    ordered_data.reverse() # reverse so order is high to low.
    ordered_data.remove(root_node) # remove the root node so it can be rsinserted at the top of the list as a source.
    ordered_data.insert(0,root_node) # place the root node at the top of the list

    path_dict = {} # to store paths
    for x in range( 0, len(ordered_data) - 1 ):
        # x will be the source node
        for y in range( x+1, len(ordered_data)):
            # y is the target node and is selected from all the nodes downstream from x
            if ordered_data[y] not in path_dict:
                # create an entry in which to store paths
                path_dict[ordered_data[y]] = []

            if nx.has_path(Z, ordered_data[x], ordered_data[y]): # test if there is a path between the compounds
                # Link together compounds
                for path in nx.all_simple_paths(Z, source = ordered_data[x], target = ordered_data[y], cutoff = depth):
                    path_dict[ordered_data[y]].append(path) # store


    paths = [] # to store output paths

    for p in path_dict:

        if len(path_dict[p]) == 0:
            pass
        else:
            scores = []
            for poss in path_dict[p]:
                line = []
                for x in poss:
                    if ">>" in x:
                        pass
                    elif x in remove_units:
                        line.append(remove_units[x])
                    else:
                        line.append(0.0)

                score = np.polyfit(np.arange(0,len(line)),line, 1)
                scores.append(score[0]/(len(poss)))

            ind = scores.index(min(scores))
            paths.append(path_dict[p][ind])

    return paths

def pathway_optmisation(G, data_container, motifs = False, root_node  = "O=C(CO)CO", depth = 20):
    '''
    General idea: find shortest path to a node from any node upstream from it
    (when nodes are ordered in decreasing amplitude).

    Parameters
    ----------
    G: NetworkX DiGraph
        Graph in which paths will be found.
    data_container: dict
        Dictionary containing node data.
    motifs: dict or False
        {"name":{"connection_point":[], "graph": nx.DiGraph}}
        Motifs are removed for the pathway search process
    root_node: str
        The compound from which the network grows.

    Returns
    -------
    paths: list of paths
    '''

    remove_units = {} # to convert keys frpom "X M" to "X"
    for d in data_container:
        new_key = d.split(" ")[0]
        remove_units[new_key] = data_container[d]

    remove_units[root_node] = 10.0 # force root_node to be highest amplitude

    Z = G.copy() # create a copy of the network
    if 'C=O' in Z.nodes:
        Z.remove_node("C=O") # remove formaldehyde so that not all pathways pass through it.

    # order the data in the dictionary by values (low to high)
    ordered_data = [x for v,x in sorted( zip( list(remove_units.values()), [*remove_units] ) ) if x != "C=O"]
    ordered_data.reverse() # reverse so order is high to low.
    ordered_data.remove(root_node) # remove the root node so it can be rsinserted at the top of the list as a source.
    ordered_data.insert(0,root_node) # place the root node at the top of the list

    path_dict = {} # to store paths
    for x in range( 0, len(ordered_data) - 1 ):
        # x will be the source node
        for y in range( x+1, len(ordered_data)):
            # y is the target node and is selected from all the nodes downstream from x
            if ordered_data[y] not in path_dict:
                # create an entry in which to store paths
                path_dict[ordered_data[y]] = []

            if nx.has_path(Z, ordered_data[x], ordered_data[y]): # test if there is a path between the compounds
                # Link together compounds
                for path in nx.all_simple_paths(Z, source = ordered_data[x], target = ordered_data[y], cutoff = depth):
                    path_dict[ordered_data[y]].append(path) # store

    paths = [] # to store output paths

    for p in path_dict:

        if len(path_dict[p]) == 0:
            pass
        else:
            scores = []
            for poss in path_dict[p]:
                line = []
                for x in poss:
                    if ">>" in x:
                        pass
                    elif x in remove_units:
                        line.append(remove_units[x])
                    else:
                        line.append(0.0)

                score = np.polyfit(np.arange(0,len(line)),line, 1)
                scores.append(score[0]/(len(poss)))

            ind = scores.index(min(scores))
            paths.append(path_dict[p][ind])

    return paths


def remove_empty_pathways(G, compounds, root_node):
    '''
    G: networkx DiGraph
    '''
    patt = Chem.MolFromSmarts('C=C')
    G2 = nx.DiGraph()
    search_G = G.copy()

    if 'C=O' in search_G.nodes:
        search_G.remove_node('C=O')
    if '[OH-]' in search_G.nodes:
        search_G.remove_node('[OH-]')
    if 'O' in search_G.nodes:
        search_G.remove_node('O')

    remove_reactions = []
    keep_reactions = []
    for ca in compounds:
        for cb in compounds:
            if ca == cb:
                pass
            elif nx.has_path(search_G, ca,cb):
                for path in nx.all_simple_paths(search_G, source = ca,target = cb,
                                            cutoff = 15):
                    tag = True
                    count = 0
                    for p in path:
                        if '>>' in p:
                            pass
                        elif Chem.MolFromSmiles(p).HasSubstructMatch(patt):
                            pass
                        elif p in compounds:
                            pass
                        elif p == root_node:
                            pass
                        else:
                            count += 1

                    if count > 0:
                        tag = False
                    if tag:
                        nx.add_path(G2,path)

    for c in compounds:
        if c not in G2.nodes:G2.add_node(c)

    if root_node not in G2.nodes:G2.add_node(root_node)

    for c in compounds:
        if c == root_node:pass
        elif nx.has_path(G2, source = root_node, target = c):pass

        elif nx.has_path(search_G, source = root_node, target = c):
            path = nx.shortest_path(search_G, root_node, target = c)
            nx.add_path(G2, path)
        elif len(G2.edges(c)) == 0:
            path_mog = []
            for c2 in compounds:
                if nx.has_path(search_G, source = c2, target = c) and not nx.has_path(G2, source = c2, target = c):
                    path = nx.shortest_path(search_G, c2, target = c)
                    if len(path) > 1:
                        path_mog.append(path)
            if len(path_mog) > 0:
                sort_paths = sorted(path_mog, key = lambda x: len(x))
                nx.add_path(G2, sort_paths[0])

    return G2

def find_paths_in_network(G,amps,root,reaction_network, search_depth = 10):
    '''
    G: networkx DiGraph
    amps: dict
    root: str
    '''
    G1 = nx.DiGraph()

    local_network_paths = network_view.paths_from_parameters_decay_score(G, amps, root_node = root, depth = search_depth)

    for x in local_network_paths:
        e_path = network_view.node_path_to_edge_path(x, reaction_network)
        G1.add_edges_from(e_path)

    for a in amps:
        if a not in [*G1.nodes]:
            G1.add_node(a)

    return G1

def connect_loose_ends(subnet, G, initial_compounds):

    paths = []
    for n in [*subnet.nodes]:
        if len(subnet.in_edges(n)) == 0 and n not in initial_compounds:
            print(n)
            sub_path = []
            min_len = 1e100
            for n2 in initial_compounds:
                if nx.has_path(G, n2, n):
                    p = nx.shortest_path(G, source = n2, target = n, method = "dijkstra")
                    if len(p) < min_len:
                        min_len = len(p)
                        sub_path = p
                else:
                    pass
            paths.append(sub_path)

    return paths

def all_shortest_paths_between_observables(G, data_container):
    # Find paths between all observed species

    observed = [x.split(" ")[0] for x in   [*data_container]]

    paths = []
    for o in observed:
        target_list = observed[:]
        target_list.remove(o)

        for t in target_list:
            if nx.has_path(G, o, t):
                path = nx.shortest_path(G, source = o, target = t)
                paths.append(path)
            else:
                pass

    return paths

def path_species_absent_from_data(G, data_dict):

    '''
    returns a list of compounds from a pathway that were not included in the
    data.

    Parameters
    ----------
    G: networkx DiGraph
        Network of the minimal reaction mechanism.
    data_dict: dictionary
        Dictionary contianing data.

    Returns
    -------

    species_list: list
        list of smiles of species which are in G, but not in data_dict
    '''

    remove_units = [] # hack to convert keys from "X M" to "X"
    for d in data_dict:
        new_key = d.split(" ")[0]
        remove_units.append(new_key)

    species_list = []
    for n in G.nodes:
        if n in remove_units:
            pass
        elif ">>" in n:
            pass
        else:
            species_list.append(n)

    return species_list
