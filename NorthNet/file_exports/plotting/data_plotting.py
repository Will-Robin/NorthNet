def plot_dataset(data, modulated_input_key = "dihydroxyacetone_concentration/ M",
                        flow_profile_key = "dihydroxyacetone_flow_profile/ uL/h"):
    '''
    For plotting a dataset- saves plots and csv files.

    Parameters
    ----------
    data: Dataset object
        Dataset object containing the data to be plotted.

    modulated_input_key: str
        Key to the concentration value for the compound which was modulated in
        the datasets conditions dictionary.

    flow_profile_key: str
        Key flow profile for the compound which was modulated in
        the datasets flow profile attribute.

    Returns
    -------
    None

    '''

    '''Figure size'''

    os.makedirs(data.name, exist_ok = True)
    os.chdir(data.name)

    up = 6
    across = 8
    # Create figure and a plotting axis
    fig = plt.figure(figsize = (12,up+5))
    ax = plt.subplot(111)

    # add second right hind y axis on which to plot the input flow profile.
    ax2 = ax.twinx()
    modulated_input = data.conditions[modulated_input_key]*data.input_flows[flow_profile_key]/data.net_flow
    ax2.plot(data.time, modulated_input*1000,"--", linewidth = 6, c = "k", alpha = 0.1)

    ordered = sorted([*data.dependents], key=lambda x: x.count("C"))
    for c in ordered:
        if c[:-2] in info_params.colour_assignments:
            ax.plot(data.time, data.dependents[c]*1000, label = info_params.smiles_to_names[c.split(" ")[0]], linewidth = 6, c  = info_params.colour_assignments[c[:-2]])
        else:
            ax.plot(data.time, data.dependents[c]*1000, label = info_params.smiles_to_names[c.split(" ")[0]], linewidth = 6)


    # Set axes parameters
    ax.tick_params(axis='both', which='major', labelsize = info_params.labels)
    ax.set_xlabel("{}".format('time/ s'), fontsize = info_params.font)
    ax.set_ylabel("concentration/ mM", fontsize = info_params.font)

    ax2.tick_params(axis='both', which='major', labelsize = info_params.labels)
    ax2.set_ylabel("{}/ mM".format(modulated_input_key.split("/")[0]), fontsize = info_params.font)
    ax2.spines['right'].set_color('grey')
    ax2.tick_params(axis='y', colors='grey')
    ax2.yaxis.label.set_color('grey')
    ax2.title.set_color('grey')

    # Place the legend in a good position.
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
                     box.width - box.x0*0.15, box.height * 0.7])

    # match the right hand y axis to the left hand axis.
    minmax = ax.get_ylim()
    ax2.set_ylim(minmax[0],minmax[1])

    # Put a legend below current axis
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4),
              fancybox=True, shadow=True, ncol=3, fontsize = info_params.legend)

    plt.savefig("{}.png".format(data.name), bbox_inches = "tight", transparent = True)
    plt.clf()

    data_out = data.time
    header = "time/ s"
    for c in data.dependents:
        data_out = np.vstack((data_out, data.dependents[c]))
        header += ", {}".format(c)

    data_out = np.vstack((data_out, modulated_input))

    header += ", {}".format(modulated_input_key)

    np.savetxt("{}.csv".format(data.name), data_out.T, header = header, delimiter = ",")
    os.chdir("..")

def plot_dataset_stack(data, modulated_input_key = "dihydroxyacetone_concentration/ M",
                        flow_profile_key = "dihydroxyacetone_flow_profile/ uL/h"):

    os.makedirs(data.name, exist_ok = True)
    os.chdir(data.name)

    font_correct = 10
    for c in data.dependents:
        fig = plt.figure(figsize = (info_params.across, info_params.up/4))
        ax = plt.subplot(111)
        if c[:-2] in info_params.colour_assignments:
            trace_name = info_params.smiles_to_names[c[:-2]]
            ax.plot(data.time, data.dependents[c]*1000, label = c.split(" ")[0], linewidth = 6, c  = info_params.colour_assignments[c[:-2]])
        else:
            trace_name = c.split(" ")[0]
            ax.plot(data.time, data.dependents[c]*1000, label = c.split(" ")[0], linewidth = 6)

        # Set axes parameters
        ax.tick_params(axis='both', which='major', labelsize = info_params.labels-font_correct)
        ax.set_yticks(np.round([np.amin(data.dependents[c]*1000), np.amax(data.dependents[c]*1000)],2))
        ax.set_xlabel("{}".format('time/ s'), fontsize = info_params.font-font_correct)
        ax.set_ylabel("concentration/ mM", fontsize = info_params.font-font_correct)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3,
                         box.width - box.x0*0.15, box.height * 0.7])

        plt.savefig("{}_{}.png".format(trace_name,data.name), bbox_inches = "tight", transparent = True)
        plt.clf()
        plt.close()

    fig = plt.figure(figsize = (info_params.across, info_params.up/4))
    ax = plt.subplot(111)
    modulated_input = data.conditions[modulated_input_key]*data.input_flows[flow_profile_key]/data.net_flow
    ax.plot(data.time, modulated_input*1000,"--", linewidth = 6, c = "k", alpha = 0.1)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
                     box.width - box.x0*0.15, box.height * 0.7])
    ax.tick_params(axis='both', which='major', labelsize = info_params.labels-font_correct)
    ax.set_ylabel("{}/ mM".format(modulated_input_key.split("/")[0]), fontsize = info_params.font-font_correct)
    ax.spines['right'].set_color('grey')
    ax.tick_params(axis='y', colors='grey')
    ax.yaxis.label.set_color('grey')
    ax.title.set_color('grey')
    plt.savefig("drive_{}.png".format(data.name), bbox_inches = "tight", transparent = True)
    plt.clf()
    plt.close()
    os.chdir("..")
    return "pfft"

def plot_normalised_dataset(data, modulated_input_key = "dihydroxyacetone_concentration/ M",
                        flow_profile_key = "dihydroxyacetone_flow_profile/ uL/h"):
    '''
    For plotting a dataset- saves plots and csv files.

    Parameters
    ----------
    data: Dataset object
        Dataset object containing the data to be plotted.

    modulated_input_key: str
        Key to the concentration value for the compound which was modulated in
        the datasets conditions dictionary.

    flow_profile_key: str
        Key flow profile for the compound which was modulated in
        the datasets flow profile attribute.

    Returns
    -------
    None

    '''

    '''Figure size'''
    up = 6
    across = 8
    # Create figure and a plotting axis
    fig = plt.figure(figsize = (12,up+5))
    ax = plt.subplot(111)

    # add second right hind y axis on which to plot the input flow profile.
    ax2 = ax.twinx()
    modulated_input = data.conditions[modulated_input_key]*data.input_flows[flow_profile_key]/data.net_flow
    ax2.plot(data.time, (modulated_input-np.average(modulated_input))/(np.amax(modulated_input)-np.average(modulated_input)),"--", linewidth = 5, c = "k", alpha = 0.5)

    for c in data.dependents:
        ax.plot(data.time, (data.dependents[c]-np.average(data.dependents[c]))/(np.amax(data.dependents[c])-np.average(data.dependents[c])), label = c, linewidth = 6)

    # Set axes parameters
    ax.tick_params(axis='both', which='major', labelsize = info_params.labels)
    ax.set_xlabel("{}".format('time/ s'), fontsize = info_params.font)
    ax.set_ylabel("concentration/ mM", fontsize = info_params.font)
    ax2.set_ylabel(modulated_input_key, fontsize = info_params.font)

    # Place the legend in a good position.
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3,
                     box.width - box.x0*0.15, box.height * 0.7])

    # match the right hand y axis to the left hand axis.
    minmax = ax.get_ylim()
    ax2.set_ylim(minmax[0],minmax[1])

    # Put a legend below current axis
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3),
              fancybox=True, shadow=True, ncol=3, fontsize = info_params.legend)

    plt.savefig("{}_normalised.png".format(data.name), bbox_inches = "tight", transparent = True)
    plt.clf()

    data_out = data.time
    header = "time/ s"
    for c in data.dependents:
        data_out = np.vstack((data_out, data.dependents[c]))
        header += ", {}".format(c)

    data_out = np.vstack((data_out, modulated_input))

    header += ", {}".format(modulated_input_key)

    np.savetxt("{}_normalised.csv".format(data.name), data_out.T, header = header, delimiter = ",")

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

def plot_species_responses(amps_dict):
    print("using species_responses_contour_plot()")
    return species_responses_contour_plot(amps_dict)

def species_responses_contour_plot(amps_dict):
    '''
    Plot output form data_to_amplitude_correlations() function.
    '''
    for d in amps_dict:
        if len(amps_dict[d]["amps"]) < 3:
            pass
        else:
            try:
                fig, ax = plt.subplots(figsize=(40,40))
                X = amps_dict[d]["drive_freq"]
                Y = [x/60 for x in amps_dict[d]["residence_time"]]
                Z = [1000*x for x in amps_dict[d]["amps"]]

                # create x-y points to be used in heatmap
                xi = np.linspace(min(X), max(X), 1000)
                yi = np.linspace(min(Y), max(Y), 1000)

                # Z is a matrix of x-y values
                zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='linear', fill_value= 0)

                # Create the contour plot
                CS = ax.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow)
                ax.scatter(X,Y, c= "k",label  = "experiment")

                cbar = plt.colorbar(mappable = CS, )
                cbar.ax.tick_params(labelsize=info_params.font)
                cbar.set_label('Amplitude', rotation=270, fontsize = info_params.font, labelpad = 20)

                ax.tick_params(labelsize = info_params.font, axis = "both")
                ax.set_xlabel("period/ min", fontsize = info_params.font+2)
                ax.set_ylabel("residence time/ min", fontsize = info_params.font+2)
                ax.set_ylim(2,8)
                ax.set_xlim(2,16)
                ax.set_title("{} amplitudes".format(d))
                plt.legend()
                plt.savefig("{}_amplitude_contour_plot.png".format(d), bbox_inches = "tight", transparent = True)

                plt.clf()
                plt.close()
            except:
                pass

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
