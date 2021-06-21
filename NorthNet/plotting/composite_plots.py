import numpy as np
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

def create_ellipse(mean, cov,
                  facecolor = 'k',edgecolor = 'k', **kwargs):

    import numpy as np
    from matplotlib.patches import Ellipse

    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    ellipse = Ellipse(mean, width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, edgecolor = edgecolor, **kwargs)

    return ellipse

def confidence_ellipse(x, y, ax = None, weights = None, n_std=3.0,
                        facecolor='none', edgecolor = 'k',**kwargs):
    """
    Adapted from matplotlib documentation.

    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    import numpy as np
    from NorthNet.plotting import composite_plots

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y, aweights = weights)

    ellipse = composite_plots.create_ellipse((0, 0), cov,
                              facecolor = facecolor,
                              edgecolor = edgecolor,
                              **kwargs)


    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
