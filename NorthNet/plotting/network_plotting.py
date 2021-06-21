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
