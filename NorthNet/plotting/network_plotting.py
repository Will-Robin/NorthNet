def plot_nodes(G, ax, color = '#000000', size = 2, alpha = 1, zorder = 1):
    '''
    Parameters
    ----------

    Returns
    -------
    None
    '''
    from NorthNet.network_manipulations.networkx_ops import coordinates

    if color == 'nodewise':
        colors = [G.nodes[n]['color'] for n in G.nodes]
    else:
        colors = color

    if size == 'nodewise':
        sizes = [G.nodes[n]['size'] for n in G.nodes]
    else:
        sizes = size


    coords = coordinates.get_network_scatter(G)

    ax.scatter(coords[0], coords[1],
               s = sizes, c = colors, alpha = alpha,
               zorder = zorder)

def draw_arrow_connectors(G,ax, color = '#000000',linew = 'edgewise', alpha = 1,
                          zorder = 0, shrink_a = 0, shrink_b = 0):
    '''
    Parameters
    ----------
    G: networkx DiGraph with position info in nodes
    ax: matplotlib axis in which to place arrows.

    Returns
    -------
    None
    '''
    from matplotlib.patches import FancyArrowPatch
    if color == 'edgewise':
        color_return = lambda x:G.edges[x]['color']
    else:
        color_return = lambda x: color

    if linew == 'edgewise':
        pass
    else:
        for e in G.edges:
            G.edges[e]['weight'] = linew

    for e in G.edges:
        begin = G.nodes[e[0]]['pos']
        end = G.nodes[e[1]]['pos']

        arrow = FancyArrowPatch(begin,end,arrowstyle='-|>', path = None,
                                connectionstyle='arc',#'Angle3'
                                zorder = zorder,
                                facecolor = color_return(e),
                                edgecolor = color_return(e),
                                linewidth = G.edges[e]['weight'],
                                mutation_scale = 5,
                                shrinkA = shrink_a,
                                shrinkB = shrink_b,
                                alpha = alpha)
        ax.add_patch(arrow)

def create_network_plot(G, ax = None, edge_color = 'k'):
    '''
    For creating a network plot

    Parameters:
    G: networkx graph
    ax: matplotlib axis object.

    Returns:
    fig,ax: matplotlib axis object.
    '''

    if ax == None:
        fig,ax = plt.subplots()

    lines = n_v.get_network_coordinates(G)

    node_pos_x = []
    node_pos_y = []
    node_c = []
    for n in G.nodes:
        node_pos_x.append(G.nodes[n]['pos'][0])
        node_pos_y.append(G.nodes[n]['pos'][1])
        if 'color' in G.nodes[n]:
            node_c.append(G.nodes[n]['color'])
        else:
            node_c.append('k')

    ax.plot(lines[0],lines[1], c = edge_color)
    ax.scatter(node_pos_x, node_pos_y, c = node_c)

    return fig, ax

def overlay_network(base_G, G, name = "data_in", fig_dims = (40,40), font = 36):

    node_sizes = {} # create a dict of node sizes
    node_cmap = {} # create a dict of node sizes
    base_edge_cmap = []

    for e in base_G.edges:
        if "color" in base_G.edges[e]:
            base_edge_cmap.append(base_G.edges[e]["color"])
        else:
            #print("No colour information for edge {}. Using black.".format(e))
            base_edge_cmap.append("#000000")

    for n in G.nodes:
        if "size" in G.nodes[n]:
            node_sizes[n] = G.nodes[n]["size"]
        else:
            node_sizes[n] = 100

        if "color" in G.nodes[n]:
            node_cmap[n] = G.nodes[n]["color"]
        else:
            node_cmap[n] = "#000000"

    if "pos" in G.nodes[[*G.nodes][0]]: # Get node coordinates from nodes
        spec_coords = {n:G.nodes[n]["pos"] for n in G.nodes}
    else:
        spec_coords = {n:(randint(1,100),randint(1,100)) for n in G.nodes}


    edge_cmap = {} # create a dict of edge colours.
    edge_widths = {} # Make dictionary of edge widths
    for e in G.edges:
        if "color" in G.edges[e]:
            edge_cmap[e] = G.edges[e]["color"]
        else:
            #print("No colour information for edge {}. Using black.".format(e))
            edge_cmap[e] = "k"

        if "width" in G.edges[e]:
            edge_widths[e] = G.edges[e]["width"]
        else:
            #print("No width information for edge {}. Using 2".format(e))
            edge_widths[e] = 2

    G2 = copy.deepcopy(G)
    remove_list = []
    for n in G2.nodes:
        if ">>" in n:
            pass
        else:
            remove_list.append(n)

    for r in remove_list:
        G2.remove_node(r)

    for e in G2.edges:
        G2.remove_edge(e)

    G2_colors = []
    G2_sizes = []
    G2_pos = []
    for n in G2.nodes:
        if "color" in G2.nodes[n]:
            G2_colors.append(G2.nodes[n]["color"])
        if "size" in G2.nodes[n]:
            G2_sizes.append(G2.nodes[n]["size"])

    if "pos" in G2.nodes[[*G2.nodes][0]]: # Get node coordinates from nodes
        spec_coords_2 = {n:G2.nodes[n]["pos"] for n in G2.nodes}
    else:
        spec_coords_2 = {n:(randint(1,100),randint(1,100)) for n in G.nodes}

    fig, ax = plt.subplots(figsize=fig_dims) # Create a figure
    nx.draw_networkx_edges(base_G, {n:base_G.nodes[n]["pos"] for n in base_G.nodes}, edge_color = base_edge_cmap)
    nx.draw(G, spec_coords, # draw the graph
            nodelist = [*G.nodes],
            node_size = [v for v in node_sizes.values()],
            node_color = [v for v in node_cmap.values()],
            font_size = font,
            #font_weight = "bold",
            font_color  = "r",
            edge_color = [v for v in edge_cmap.values()],
            width = [v for v in edge_widths.values()],
            node_shape = "o", alpha = 1.0,arrowsize=30, arrowstyle='-|>')

    nx.draw(G2, spec_coords_2, # draw the graph
            nodelist = [*G2.nodes],
            node_size = G2_sizes,
            node_color = "#000000",
            font_size = font,
            font_color  = "r",
            font_weight = "bold",
            node_shape = "o")

    plt.savefig("{}_network.png".format(name), bbox_inches = "tight", transparent = True)
    plt.clf()
    plt.close()

    return None
