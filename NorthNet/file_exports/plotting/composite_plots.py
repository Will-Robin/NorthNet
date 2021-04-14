def chord_diagram(corr_mat, labelling = [], colour_map = [], dat_labels = [], ax = None):

    from NorthNet import data_processing as d_p

    if ax == None:
        fig, ax = plt.subplots()
    if len(labelling) == 0:
        labelling = {1:[x for x in range(0,len(corr_mat))]}
    if colour_map == []:
        import matplotlib
        col_map = matplotlib.cm.get_cmap('tab20c', len(corr_mat))
        colour_map = [ matplotlib.colors.rgb2hex( col_map(i)[:3] ) for i in range(col_map.N) ]
    if len(dat_labels) != len(corr_mat):
        dat_labels = ['' for x in range(0,len(corr_mat))]

    circ_rad = 1
    ang_unit = 360/len(corr_mat)
    start_angle = -ang_unit/2
    # build positioning by cluster
    clusters = [*labelling]
    angles = {}
    pie_colours = []
    re_label = []
    ang_counter = 0
    for x in range(0,len(clusters)):
        inds = labelling[clusters[x]]
        for i in inds:
            angles[i] = start_angle + ang_unit*(ang_counter+0.5)
            pie_colours.append(colour_map[x])
            re_label.append(dat_labels[i])
            ang_counter += 1

    ax.set_aspect('equal')
    visited = []
    for x in range(0,len(corr_mat)):
        for y in range(0,len(corr_mat)):
            for l in labelling:
                if x in labelling[l]:
                    clust_ind = l
                    break
            if y not in labelling[clust_ind] and y not in visited and corr_mat[x,y] > 0.2:

                x1,y1 = d_p.polar_to_xy(circ_rad, angles[x])
                x2,y2 = d_p.polar_to_xy(circ_rad, angles[y])

                source = [x1,y1]
                target = [x2,y2]
                centre = [0,0]
                xnew,ynew = d_p.curve_between_points(source, target, centre)

                ax.plot(xnew, ynew, linewidth = abs(corr_mat[x,y])/2, zorder = 0,
                        c = colour_map[clust_ind-1], alpha = abs(corr_mat[x,y]))

                visited.append(x)

    wedges,texts = ax.pie([1 for v in range(0,len(corr_mat))],
                    explode = [1/len(corr_mat) for x in range(0,len(corr_mat))],
                    colors = pie_colours, startangle = start_angle,
                    labels = re_label,
                    radius = circ_rad,
                    rotatelabels =True,
                    textprops = {'fontsize':3})

    for p in wedges:
        p.set_linewidth(0.2)
        p.set_edgecolor('k')
        p.set_width(1/10)

    return ax

def chord_diagram_2(corr_mat, labelling = [], colour_map = [], dat_labels = [],
                    ax = None, radius = 1, fig = None):

    from NorthNet import data_processing as d_p

    if ax == None:
        ax = plt.subplot(111, projection='polar')
    if len(labelling) == 0:
        labelling = {1:[x for x in range(0,len(corr_mat))]}
    if colour_map == []:
        import matplotlib
        col_map = matplotlib.cm.get_cmap('tab20c', len(corr_mat))
        colour_map = [ matplotlib.colors.rgb2hex( col_map(i)[:3] ) for i in range(col_map.N) ]
    if len(dat_labels) != len(corr_mat):
        dat_labels = ['' for x in range(0,len(corr_mat))]

    inner_radius = radius*0.9
    radii = [radius-inner_radius for x in range(0,len(corr_mat))]
    width = 2 * np.pi/len(corr_mat)
    theta = np.zeros(len(corr_mat))
    colours = [0]*len(corr_mat)
    t_count = 0
    for l in labelling:
        for x in labelling[l]:
            theta[x] = t_count*(2*np.pi)/len(corr_mat)
            colours[x] = colour_map[l-1]
            t_count += 1
    visited = []

    for l in labelling:
        clust_ind = l
        for x in labelling[l]:
            for y in range(0,len(corr_mat)):
                if y not in labelling[clust_ind] and y not in visited and corr_mat[x,y] > 0.9:

                    x_new = np.linspace(theta[x],theta[y],num=30)
                    y_new = inner_radius*np.cos(x_new+theta[x])

                    ax.plot(x_new,y_new,
                            linewidth = abs(corr_mat[x,y])/2, zorder = 2,
                            c = colour_map[clust_ind-1], alpha = abs(corr_mat[x,y]))

                    visited.append(x)

    bars = ax.bar(theta, radii, width=width, bottom= inner_radius,
                  color = colours, alpha = 1, edgewidths = 0)

    return ax

def draw_radar(x_axis, y_axis,base = 0,ax = None, colour = 'k'):
    if ax == None:
        ax = plt.subplot(111, projection='polar')

    width = 2 * np.pi/len(x_axis)
    width = width + width/100
    theta = np.linspace(0.0, 2 * np.pi, len(x_axis), endpoint=False)
    radii = y_axis

    bars = ax.bar(theta, radii, width=width, bottom=base, edgecolor = 'w',
                  color = colour, alpha = 1, linewidth = 0.0)

    return ax

def network_in_expr_circle(circ_data,circ_col_map, pie_rad, network, reactotypes,
                            network_linew,  network_coords, base_network = [],
                            base_net_alpha = 0.6, ax = None,
                            chord_colour = '#000000'):

    from NorthNet import conversions as conv
    from NorthNet import network_view as n_v
    n_removals = ["C=O","O","[OH-]"]

    if ax == None:
        fig,ax = plt.subplots()
    if len(circ_col_map) != len(circ_data):
        circ_col_map = ['w' for x in circ_data]

    wedges,texts = ax.pie(circ_data, colors = circ_col_map, radius = pie_rad)
    for p in wedges:
        p.set_linewidth(0.25)
        p.set_edgecolor('k')
        p.set_width(pie_rad/3)

    srt = sorted(network, key = lambda x:circ_data[[*reactotypes].index(x.Name)])
    if len(base_network) != 0:
        ax.plot(base_network[0], base_network[1], '--',c= '#000000', zorder = 0, alpha = base_net_alpha, linewidth = network_linew/2)

    if len(network) != 0:
        for c2,col in enumerate(srt):

            if len(col.NetworkReactions) == 0.0:
                continue

            G = conv.convert_to_networkx(col)

            for node in n_removals:
                if node in G.nodes:
                    G.remove_node(node)

            for node in G.nodes:
                G.nodes[node]['pos'] = network_coords[node]

            xy = n_v.get_network_coordinates(G)

            alph = circ_data[[*reactotypes].index(col.Name)]

            ax.plot(xy[0], xy[1], c=chord_colour, zorder = c2+1, alpha = alph, linewidth = network_linew)

    ax.set_axis_off()
    ax.set_aspect('equal')

    return ax
