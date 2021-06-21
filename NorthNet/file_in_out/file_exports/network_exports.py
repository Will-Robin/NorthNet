def save_networkx_graph(G, name = "saved_networkx_graph"):
    '''
    '''
    Z = G.copy()
    for n in Z.nodes:
        if "pos" in Z.nodes[n]:
            Z.nodes[n]["pos"] = "coordinate removed for saving"
        if ">>" in n:
            Z.nodes[n]["type"] = "reaction"
        else:
            Z.nodes[n]["type"] = "compound"

    nx.write_gexf(Z, "{}.gexf".format(name))

def pickle_network_data(network, name):
    '''
    For pickling network data. Effectively a wrapper for pickle.
    Parameters
    ----------
    network: NorthNet_functions Network object (or any other object)
        Object to be pickled.
    name: str
        A name for the output file.

    Returns
    -------
    None
    '''

    with open("{}_pickled".format(name), "wb") as f:
        pickle.dump(network,f)

def network_to_spreadsheet(network,fname):
    '''
    Places the reactions in a network into a spreadsheet (comma separated values).
    Parameters
    ----------
    network: NorthNet Network object
        Network to be converted to spreadheet.
    fname: str
        A name for the output file.
    Returns
    -------
    None
    '''
    with open('{}.txt'.format(fname), 'w') as f:
        [f.write('{}\t'.format(h)) for h in network.NetworkReactions[ [*network.NetworkReactions][0] ][0].Database_Entries[0].get_header() ]
        f.write('\n')
        for reaction in [*network.NetworkReactions]:
            for x in network.NetworkReactions[reaction].Database_Entries:
                [f.write('{}\t'.format(y)) for y in x.recompile()]
                f.write('\n')

def write_as_networkx(network, extension = "gexf", sv_img = False):
    '''
    Creates a reaction-species graph from a network object using networkx.
    The extension is the filetype to be output.

    Parameters
    ----------
    network: NorthNet_functions Network object
        NorthNet Network to be converted to networkx network.
    extension: str
        extension for the output file.
    sv_img: Bool
        Whether to create images in the file output or not.

    Returns
    -------
    fname: str
        Name of the output file.
    '''

    fname = "{}.{}".format(network.Name,extension)
    G = conv.convert_to_networkx(network, save_images = sv_img)
    nx.write_gexf(G, fname)

    return fname

def SNetwork_write_as_networkx(network, extension = "gexf", sv_img = False):
    '''
    Creates a reaction-substructure-species graph from a network object using networkx.
    The extension is the filetype to be output.

    Parameters
    ----------
    network: NorthNet_functions SNetwork object
        NorthNet SNetwork to be converted to networkx network.
    extension: str
        extension for the output file.
    sv_img: Bool
        Whether to create images in the file output or not.

    Returns
    -------
    fname: str
        Name of the output file.
    '''

    fname = "{}.{}".format(network.Name,extension)
    G = conv.SNetwork_convert_to_networkx(network, save_images = sv_img)
    nx.write_gexf(G, fname)

    return fname

def graphviz_output_old(network, fname, render_engine = 'fdp'):
    '''
    Creates a graphviz graph from a NorthNet Network.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be represented.
    fname: str
        A name for the output file.
    render_engine: str
        Layout render engine for graphviz.
    Returns
    -------
    None
    '''

    dot = Digraph(comment = fname, engine=render_engine, strict = 'True') # Create a graph with graphviz to plot a scheme of the network
    dot.attr(newrank = 'true')
    #dot.attr(rankdir = 'LR',size='1000,1000') # removed if complications arise in graph creation. Otherwise,m it should help line up species

    picturenames = [] # a list of names for generated picture files

    cont = 0
    for n in network.NetworkCompounds:
        if Chem.MolFromSmiles(n).GetNumAtoms() == 1:
            pass
        else:
            A_name = drw.create_image(n, str(cont))
            picturenames.append(A_name)
            dot.node(n, "",image = '{}.png'.format(cont),
            style='filled',
            shape = 'circle',
            color = "blue")
            cont += 1

    for r in network.NetworkReactions:    # Create edges between the nodes.

        reactants = network.NetworkReactions[r].Reactants
        products = network.NetworkReactions[r].Products
        '''
        im = drw.draw_reaction(network.NetworkReactions[r].ReactionSMILES,network.NetworkReactions[r].ReactionSMILES)
        picturenames.append(im)
        dot.node(network.NetworkReactions[r].ReactionSMILES, "",image = '{}'.format(im),
        style='filled',
        shape = 'square',
        color = "grey")
        '''
        r_node_name = network.NetworkReactions[r].Database_Entries[0].Info["Reaction ID"]
        dot.node(r_node_name, r_node_name)

        for reac in network.NetworkReactions[r].Reactants:
            if network.NetworkCompounds[reac].Mol.GetNumAtoms() == 1:
                pass
            else:
                dot.edge(reac,r_node_name)

        for p in network.NetworkReactions[r].Products:
            if network.NetworkCompounds[p].Mol.GetNumAtoms() == 1:
                pass
            else:
                dot.edge(r_node_name,p)

    dot.render("{}.gv".format(fname), view = False, renderer = "cairo") # render the graph generated

    for f in picturenames: # delete the possibly large amount of images generated for the nodes
        os.remove(f)

    return 0

def graphviz_output_substructure(network, fname, render_engine = 'fdp'):
    '''
    Creates a graphviz graph from a NorthNet SNetwork.

    Parameters
    ----------
    network: NorthNet SNetwork object
        NorthNet SNetwork to be represented.
    fname: str
        A name for the output file.
    render_engine: str
        Layout render engine for graphviz.
    Returns
    -------
    None
    '''

    dot = Digraph(comment = fname, engine=render_engine, strict = 'True') # Create a graph with graphviz to plot a scheme of the network
    dot.attr(newrank = 'true')

    picturenames = [] # a list of names for generated picture files

    cont = 0
    for n in network.SNetworkSubstructs:
        c_node_name  = "S_"+Chem.MolToInchi(network.SNetworkSubstructs[n].Mol)

        A_name = drw.create_image(n, str(cont), fig_size = (150,150))
        picturenames.append(A_name)
        dot.node(c_node_name, "",image = '{}.png'.format(cont),
        style='filled',
        shape = 'circle',
        color = "red")
        cont += 1

    for n in network.SNetworkCompounds:
        c_node_name  = "C_"+Chem.MolToInchi(network.SNetworkCompounds[n].Mol)

        A_name = drw.create_image(n, str(cont), fig_size = (150,150))
        picturenames.append(A_name)
        dot.node(c_node_name, "",image = '{}.png'.format(cont),
        style='filled',
        shape = 'circle',
        color = "blue")
        cont += 1

        for i in network.SNetworkCompounds[n].In:
            dot.edge( "S_"+ Chem.MolToInchi(network.SNetworkSubstructs[i].Mol), c_node_name )
        for o in network.SNetworkCompounds[n].Out:
            dot.edge(c_node_name, "S_"+ Chem.MolToInchi(network.SNetworkSubstructs[o].Mol) )

    for r in network.SNetworkTemplates: # Create edges between reactions and substructures.
        if len(network.SNetworkTemplates[r].Database_Entries) == 0:
            r_node_name = network.SNetworkTemplates[r].Generation_Details[0].ReactionSMILES
        else:
            r_node_name = network.SNetworkTemplates[r].Database_Entries[0].Info["Reaction ID"]

        dot.node(r_node_name, r_node_name)

        for reac in network.SNetworkTemplates[r].ReactionTemplate.ReactantSubstructures:
            dot.edge("S_"+ Chem.MolToInchi(network.SNetworkSubstructs[reac].Mol), r_node_name)

        for p in network.SNetworkTemplates[r].ReactionTemplate.ProductSubstructures:
            dot.edge(r_node_name, "S_"+Chem.MolToInchi(network.SNetworkSubstructs[p].Mol))

    dot.render("{}_SNetwork.gv".format(fname), view = False, renderer = "cairo") # render the graph generated

    for f in picturenames: # delete the possibly large amount of images generated for the nodes
        os.remove(f)

    return 0

def write_reaction_rules_from_network(network,fname):
    '''
    network: NorthNet Network Object
    '''

    r_list = []
    for r in network.NetworkReactions:
        r_list.append(network.NetworkReactions[r])

    substruct_net = Classes.Substructure_Network(r_list, "{}_substruct".format(''), "")

    # write new reaction rules to .csv file
    with open(fname,'w') as f:
        f.write('Reaction Name\tReactant_substructures\tProduct_substructures\tReaction SMARTS\n')
        for r in substruct_net.SNetworkTemplates:
            r_temp = substruct_net.SNetworkTemplates[r]
            f.write('{}\t{}\t{}\t{}'.format(r_temp.ReactionTemplate.Name, '.'.join(r_temp.ReactionTemplate.ReactantSubstructures), '.'.join(r_temp.ReactionTemplate.ProductSubstructures), r_temp.ReactionTemplate.ReactionSMARTS))
            f.write('\n')

def graph_by_correlation(groups, amplitudes, graph_name):

    G = nx.DiGraph()
    for g in groups:
        specs = g.split(",")
        if specs[0] == specs[1]:
            pass
        else:
            G.add_node(specs[0], label = specs[0], size_val = float(amplitudes[specs[0]]))
            G.add_node(specs[1], label = specs[1], size_val = float(amplitudes[specs[1]]))
            if float(groups[g]) > 0.7:
                G.add_edge(specs[1],specs[0], Weight  = float(groups[g]))
            elif float(groups[g]) < -0.7:
                G.add_edge(specs[0],specs[1], Weight  = float(groups[g]))
            else:
                pass

    nx.write_gexf(G, "{}.gexf".format(graph_name))
    print('Made graph.')