def plot_nodes(G, ax, color = '#000000', size = 2, alpha = 1, zorder = 1):
    '''
    Parameters
    ----------

    Returns
    -------
    None
    '''
    from NorthNet.networkx_manipulations.networkx_ops import coordinates

    if color == 'nodewise':
        colors = [G.nodes[n]['color'] for n in G.nodes]
    else:
        colors = color

    if size == 'nodewise':
        sizes = [G.nodes[n]['size'] for n in G.nodes]
    else:
        sizes = size


    coords = coordinates.get_network_scatter(G)
    print(coords.shape)
    print(len(colors))

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

def plot_network(G, node_labels = [], name = "data_in", second_node_param_set = [],
                fig_dims = (40,40), font = 36, background_color = '#000000',
                annotations = []):
    '''
    Add data into networkx object

    Parameters
    ----------
    G: networkx.DiGraph object
        Data will be inserted into this Networkx DiGraph.

    name: str
        Name for the output graph picture.

    Returns
    -------
    G: networkx.DiGraph object
        Original graph.
    '''


    if len(node_labels) == 0: # Create dictionary for labels
        lbs = {}#lbs = {n:n for n in G.nodes} # if no specific labels are provided, label everything
    else:
        remove_units = [x.split(" ")[0] for x in node_labels]
        blank_names = [x.split("_")[0] for x in remove_units]
        lbs = {s:b for b,s in zip(blank_names,remove_units)} # Create labels for graph from node_labels.

    node_sizes = {} # create a dict of node sizes
    node_cmap = {} # create a dict of node sizes
    for n in G.nodes:
        if "size" in G.nodes[n]:
            node_sizes[n] = G.nodes[n]["size"]
        else:
            node_sizes[n] = 100

        if "color" in G.nodes[n]:
            node_cmap[n] = G.nodes[n]["color"]
        else:
            node_cmap[n] = "k"

    if "pos" in G.nodes[[*G.nodes][0]]: # Get node coordinates from nodes
        spec_coords = {n:G.nodes[n]["pos"] for n in G.nodes}
    else:
        spec_coords = {n:(randint(1,100),randint(1,100)) for n in G.nodes}

    remove_background_edges = []
    base_G = nx.DiGraph()
    for e in G.edges:
        if G.edges[e]["color"] == background_color:
            base_G.add_edge(e[0],e[1], color = background_color, width = G.edges[e]["width"])
            remove_background_edges.append(e)

    for e in remove_background_edges:
        G.remove_edge(e[0],e[1])

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
    nx.draw_networkx_edges(base_G, {n:G.nodes[n]["pos"] for n in base_G.nodes}, edge_color = background_color, width =10.0, alpha = 0.5)
    nx.draw(G, spec_coords, # draw the graph
            nodelist = [*G.nodes],
            node_size = [v for v in node_sizes.values()],
            node_color = [v for v in node_cmap.values()],
            font_size = font,
            #font_weight = "bold",
            font_color  = "r",
            labels = lbs,
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

    if len(second_node_param_set) > 0:
        for n in second_node_param_set:
            if n in G.nodes:
                ax.scatter(G.nodes[n]["pos"][0],G.nodes[n]["pos"][1], facecolors="None", edgecolors = G.nodes[n]["color"],  s =second_node_param_set[n], linewidth = 7, zorder = 1)

    if len(annotations) > 0:
        for a in annotations:
            plt.annotate(a, xy = annotations[a], fontsize = font)

    plt.savefig("{}_network.png".format(name), bbox_inches = "tight", transparent = True)
    plt.clf()
    plt.close()

    return G, G2

def plot_network_2(G):
    import plotly.graph_objects as go
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=1 , color='#888'),
                hoverinfo='none',
                mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    G_colors = []
    G_sizes = []
    for n in G.nodes:
        if "color" in G.nodes[n]:
            G_colors.append(G.nodes[n]["color"])
        if "size" in G.nodes[n]:
            G_sizes.append(G.nodes[n]["size"])

    node_trace = go.Scatter(
                            x=node_x, y=node_y,
                            mode='markers',
                            hoverinfo='text',
                            marker=dict(
                                showscale=True,
                                # colorscale options
                                #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
                                #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
                                #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                                colorscale='YlGnBu',
                                reversescale=True,
                                color=G_colors,
                                size=G_sizes,

                                colorbar=dict(
                                    thickness=15,
                                    title='Node Connections',
                                    xanchor='left',
                                    titleside='right'
                                ),

                                line_width=2))

    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='<br>Network graph made with Python',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    annotations=[ dict(
                        text="Python code: <a href='https://plot.ly/ipython-notebooks/network-graphs/'> https://plot.ly/ipython-notebooks/network-graphs/</a>",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.show()

def put_images_in_nodes(G, node_labels = [], name = "data_in"):
    '''
    Add data into networkx object

    Parameters
    ----------
    G: networkx.DiGraph object
        Data will be inserted into this Networkx DiGraph.

    name: str
        Name for the output graph picture.

    Returns
    -------
    G: networkx.DiGraph object
        Original graph.
    '''


    if len(node_labels) == 0: # Create dictionary for labels
        lbs = {} # if no specific labels are provided, label everything
    else:
        remove_units = [x.split(" ")[0] for x in node_labels]
        lbs = {s:s for s in remove_units} # Create labels for graph from node_labels.

    node_sizes = {} # create a dict of node sizes
    for n in G.nodes:
        if "size" in G.nodes[n]:
            if ">>" in n:
                node_sizes[n] = 0
            else:
                node_sizes[n] = G.nodes[n]["size"]

        else:
            #print("No size information for node {}.".format(n))
            #print("Using size 10.")
            node_sizes[n] = 10

    node_cmap = {} # create a dict of node sizes
    for n in G.nodes:
        if "color" in G.nodes[n]:
            node_cmap[n] = G.nodes[n]["color"]
        else:
            #print("No color information for node {}.".format(n))
            #print("Using black.")
            node_cmap[n] = "k"

    edge_cmap = {} # create a dict of edge colours.
    for e in G.edges:
        if "color" in G.edges[e]:
            edge_cmap[e] = G.edges[e]["color"]
        else:
            #print("No colour information for edge {}. Using black.".format(e))
            edge_cmap[e] = "k"

    edge_widths = {} # Make dictionary of edge widths
    for e in G.edges:
        if "width" in G.edges[e]:
            edge_widths[e] = G.edges[e]["width"]
        else:
            #print("No width information for edge {}. Using 2".format(e))
            edge_widths[e] = 2

    if "pos" in G.nodes[[*G.nodes][0]]: # Get node coordinates from nodes
        spec_coords = {n:G.nodes[n]["pos"] for n in G.nodes}
    else:
        #print("No coordinate information in nodes. Key: ""pos"".")
        #print("Using random coordinates instead.")
        spec_coords = {n:(randint(1,100),randint(1,100)) for n in G.nodes}

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
        #print("No coordinate information in nodes. Key: ""pos"".")
        #print("Using random coordinates instead.")
        spec_coords_2 = {n:(randint(1,100),randint(1,100)) for n in G.nodes}
    '''
    Add images into graph nodes
    '''
    file_list = os.listdir()
    for n in G.nodes:

        if "{}.png".format(n) in file_list:
            G.nodes[n]["image"] = "{}.png".format(n)
        elif "S:" in n or ">>" in n:
            G.nodes[n]["image"] = "None"
        else:
            fname = drawing.create_image(n, n, fig_size = (300,300))
            G.nodes[n]["image"] = fname

    fig, ax = plt.subplots(figsize=(40,40)) # Create a figure

    for n in G:
        if ">>" in n or "S:" in n:
            pass
        elif G.nodes[n]['image'] != "None":
            xx,yy = G.nodes[n]["pos"]
            img = Image.open(G.nodes[n]['image'])
            x = np.asarray(img.convert('RGBA')).copy()

            x[:, :, 3] = (255 * (x[:, :, :3] == 0).any(axis=2)).astype(np.uint8)
            newimg = Image.fromarray(x)

            newimg = newimg.resize((int(0.5*len(x)),int(0.5*len(x))))
            ab = AnnotationBbox(OffsetImage(newimg), (xx,yy), frameon=False)
            ax.add_artist(ab)
        else:
            pass

    nx.draw(G, spec_coords, # draw the graph
            nodelist = [*G.nodes],
            node_size = [v for v in node_sizes.values()],
            node_color = [v for v in node_cmap.values()],
            font_size = 16,
            font_weight = "bold",
            labels = lbs,
            edge_color = [v for v in edge_cmap.values()],
            width = [v for v in edge_widths.values()],
            node_shape = "o", alpha = 0.4,arrowsize=30, arrowstyle='-|>')

    nx.draw(G2, spec_coords_2, # draw the graph
            nodelist = [*G2.nodes],
            node_size = G2_sizes,
            node_color = G2_colors,
            font_size = 16,
            font_weight = "bold",
            node_shape = "d")

    plt.savefig("{}_network.png".format(name), transparent = True)
    plt.clf()
    plt.close()

    return G
