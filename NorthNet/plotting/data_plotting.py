def create_ellipse(mean, cov,
                  facecolor = 'k',edgecolor = 'k', **kwargs):

    from matplotlib.patches import Ellipse
    import numpy as np

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
    from NorthNet.file_exports.plotting import data_plotting

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y, aweights = weights)

    ellipse = data_plotting.create_ellipse((0, 0), cov,
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

def clustering_diagram_3D(Y,labels, groups,graph_name):
    '''
    Create a clustering diagram in 3D from correlations between species.

    Y:  2D numpy array
        Euclidean distance matrix
    groups: dict
        dictionary indexed by variables paired by correlation values.
    graph_name: str
        a name for the output plot.


    '''

    fig = plt.figure(figsize=(20,20))
    ax = fig.add_subplot(111, projection='3d')
    all_corrs = [abs(x) for x in groups.values()]

    for g in groups:
        line = [int(x) for x in g.split(",")]
        ax.scatter([ Y[line[0],0], Y[line[1],0] ],[ Y[line[0],1], Y[line[1],1] ], [ Y[line[0],2], Y[line[1],2] ], marker = ".", s = 300)
        ax.text(Y[line[0],0], Y[line[0],1], Y[line[0],2], labels[line[0]], fontsize = info_params.font-5 , ha = "left")
        ax.text(Y[line[1],0], Y[line[1],1], Y[line[1],2], labels[line[1]], fontsize = info_params.font-5 , ha = "left")

        if abs(groups[g]) > 0.5:
            ax.plot([ Y[line[0],0], Y[line[1],0] ],[ Y[line[0],1], Y[line[1],1] ], [ Y[line[0],2], Y[line[1],2] ],linewidth = abs(groups[g])*5)



    ax.tick_params(labelsize = info_params.font, axis = "both")
    plt.savefig("{}_clustering_3D.png".format(graph_name), bbox_inches = "tight", transparent = True)
    plt.clf()
    plt.close()

def clustering_diagram_2D(Y,labels, groups,graph_name):
    '''
    Create a clustering diagram in 2D from correlations between species.
    Y:  2D numpy array
        Euclidean distance matrix
    groups: dict
        dictionary indexed by variables paired by correlation values.
    graph_name: str
        a name for the output plot.

    '''


    all_corrs = [abs(x) for x in groups.values()]

    visited = []
    for pc1 in range(0,len(Y[0])):
        for pc2 in range(0,len(Y[0])):
            if pc1 in visited or pc1 == pc2:
                pass
            else:
                fig = plt.figure(figsize=(info_params.across,info_params.up))
                ax = fig.add_subplot(111)
                for g in groups:
                    line = [labels.index(x) for x in g.split(",")]
                    c1 = [x.split(" ")[0] for x in g.split(",")]
                    colors = [info_params.colour_assignments[x] for x in c1]
                    ax.scatter([ Y[line[0],pc1], Y[line[1],pc1] ],[ Y[line[0],pc2], Y[line[1],pc2] ], marker = ".", s = 500, linewidth = groups[g], c = colors)
                    #ax.text(Y[line[0],0], Y[line[0],1], labels[line[0]], fontsize = info_params.font-20 , ha = "left", va = "top")
                    #ax.text(Y[line[1],0], Y[line[1],1], labels[line[1]], fontsize = info_params.font-15 , ha = "center", va = "top")

                    if abs(groups[g]) > 0.8:
                        ax.plot([ Y[line[0],pc1], Y[line[1],pc1] ],[ Y[line[0],pc2], Y[line[1],pc2] ], alpha = 0.5, c = "k", zorder = 0)

                ax.plot([np.amin(Y),np.amax(Y)],[0,0],"--",c = "k", alpha = 0.6)
                ax.plot([0,0],[np.amin(Y),np.amax(Y)],"--",c = "k", alpha = 0.6)
                ax.tick_params(labelsize = info_params.labels, axis = "both")
                ax.set_xlabel("PC{}".format(pc1),fontsize = info_params.font)
                ax.set_ylabel("PC{}".format(pc2),fontsize = info_params.font)
                box = ax.get_position()
                ax.set_position([box.x0+ box.x0*0.6, box.y0 + box.height * 0.1,
                                 box.width - box.x0*0.6, box.height * 0.9])
                plt.savefig("{}_PC_{}_{}_clustering_2D.png".format(graph_name,pc1+1,pc2+1), bbox_inches = "tight", transparent = True)
                plt.clf()
                plt.close()

                visited.append(pc1)

def colorgrid(matrix,deps,cbar_label, cmap_name, graph_name):

    fig, ax = plt.subplots(figsize=(25,25))
    ax.set_xticklabels(deps, fontsize = info_params.font, rotation = 45, ha = "right", va = "top", color= "k")
    ax.set_yticklabels(deps, fontsize = info_params.font, rotation = 45, ha = "left",va = "center",color= "k", horizontalalignment = "center")
    ax.set_xticks(np.arange(0,len(deps)))
    ax.set_yticks(np.arange(0,len(deps)))
    ax.tick_params(labelsize=info_params.font, direction = "in", length = 5, width = 5)
    plt.imshow(matrix, cmap = cmap_name)
    try:
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=info_params.font, direction = "in", length = 5, width = 5)
        cbar.set_label(cbar_label, rotation=270, fontsize = info_params.font, labelpad = 20)
    except:
        print("Colourbar issue, colorgrid function, plotting_operations")
    plt.savefig("{}.png".format(graph_name), bbox_inches = "tight", transparent = True)
    plt.clf()
    plt.close()

def plot_cycle(amplitude_dict, time_lag_dict, name = "name"):

    '''
    Parameters
    -------------
    amplitude_dict: dict

    time_lag_dict: dict

    '''
    fig, ax = plt.subplots(figsize=(10,10))
    max_H = 0
    for a in amplitude_dict:
        H = amplitude_dict[a]
        phi = time_lag_dict[a]

        X = H*np.sin(phi)
        Y = H*np.cos(phi)

        ax.plot([0,X],[0,Y], linewidth = info_params.lines)
        ax.scatter(X,Y)

        ax.annotate(a[:-2], xy = (X,Y))

        if H > max_H:
            max_H = H

    circle2 = plt.Circle( (0,0), radius = max_H, color = "k", fill=False)
    ax.add_artist(circle2)
    ax.set_xlim(-1.1*max_H, 1.1*max_H)
    ax.set_ylim(-1.1*max_H, 1.1*max_H)
    plt.savefig("{}.png".format(name), bbox_inches = "tight", transparent = True)
    plt.clf()

def plot_trajectories(d_sets, period_sequence):
    '''

    '''
    x_vals = {}
    y_vals = {}
    for d in d_sets:
        for dep in d_sets[d].dependents:
            if dep in x_vals:
                pass
            else:
                x_vals[dep] = []
                y_vals[dep] = []

    max_H = 0
    for p,d in enumerate(d_sets):
        os.makedirs("{}".format(d_sets[d].name), exist_ok = True)
        os.chdir("{}".format(d_sets[d].name))

        tlcm = d_p.time_lag_corr_mat(d_sets[d])

        time_lags = ((d_p.time_lags(tlcm,d_sets[d].time)*360/(period_sequence[p]*60))%360)

        tl_to_dr = ((d_p.time_lags_to_drive(tlcm,d_sets[d].time)*360/(period_sequence[p]*60))%360)
        time_lag_dict = {d:t for d,t in zip([*d_sets[d].dependents], tl_to_dr)}
        amps = d_p.get_amplitudes(d_sets[d], write_f = False)


        for a in amps:
            H = amps[a]
            phi = time_lag_dict[a]
            if H > max_H:
                max_H = H

            x_vals[a].append(H*np.sin(phi))
            y_vals[a].append(H*np.cos(phi))

    fig, ax = plt.subplots(figsize=(10,10))
    for x in x_vals:
        ax.plot(x_vals[x],y_vals[x], linewidth = info_params.lines)
        ax.scatter(x_vals[x][-1],y_vals[x][-1])

        ax.annotate(x[:-2], xy = (x_vals[x][-1],y_vals[x][-1]))

    circle2 = plt.Circle( (0,0), radius = max_H, color = "k", fill=False)
    ax.scatter(0,0, c = "k")
    ax.add_artist(circle2)
    ax.set_xlim(-1.1*max_H, 1.1*max_H)
    ax.set_ylim(-1.1*max_H, 1.1*max_H)
    plt.show()
