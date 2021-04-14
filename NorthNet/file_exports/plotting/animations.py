def generate_timecourse_gif(dataset, G, out_name, flow_key = "DHA_flow/ uL/h",
                                pass_fraction = 10, scale_factor = 1e5):
    '''
    Parameters
    -----------

    dataset: NetFit Dataset object
        Dataset for plotting
    G: networkx Digraph object
        Graph for representation
    out_name: str
        Name for output gif.
    flow_key: str
        key to flow profile to be displayed in plot inset (in dataset.input_flows)
    pass_fraction: int
        Fraction of frames to take from timecourse
    scale_factor: float
        Factor by which to scale the data for plotting.

    Returns
    -------
    None
    '''
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from imageio import mimsave
    
    images = [] # array to store images in
    image_names = []

    data_container = {}
    for d in dataset.dependents:
        ins = d.split(" ")[0]
        data_container[ins] = dataset.dependents[d]

    lbs = {s:s for s in [*data_container]}

    # Get node coordinates from nodes
    if "pos" in G.nodes[[*G.nodes][0]]:
        spec_coords = {n:G.nodes[n]["pos"] for n in G.nodes}
    else:
        print("No coordinate information in nodes. Key: ""pos"".")
        print("Using random coordinates instead.")
        spec_coords = {n:(randint(1,100),randint(1,100)) for n in G.nodes}

    # create a list of node colours.
    cmap = []
    for n in G.nodes:
        if "color" in G.nodes[n]:
            cmap.append(G.nodes[n]["color"])
        else:
            print("No colour information for node {}. Using black.".format(n))
            cmap.append("k")

    # create a list of edge colours.
    cmap2 = []
    for e in G.edges:
        if "color" in G.edges[e]:
            cmap2.append(G.edges[e]["color"])
        else:
            print("No colour information for edge {}. Using black.".format(e))
            cmap2.append("k")


    for t in range(0,len(dataset.time)):

        if t % pass_fraction == 0:

            for s in G.nodes:
                if s in data_container:
                    G.nodes[s]["size"] = data_container[s][t]*scale_factor
                else:
                    G.nodes[s]["size"] = 1

            # Get node sizes from nodes.
            if "size" in G.nodes[[*G.nodes][0]]:
                node_sizes = {n:G.nodes[n]["size"] for n in G.nodes}
            else:
                print("No size information in nodes. Key: ""size"".")
                print("Using size 10.")
                node_sizes = {n:10 for n in G.nodes}

            # Create a figure
            fig = plt.figure(figsize = (20,20))
            ax1 = fig.add_subplot(111)

            # draw the graph
            nx.draw(G, spec_coords,
                    nodelist = [*G.nodes],
                    node_size = [v for v in node_sizes.values()],
                    node_color = cmap,
                    font_size = 16,
                    font_weight = "bold",
                    labels = lbs,
                    edge_color = cmap2,
                    width = 2)

            axins = inset_axes(ax1, width=5, height=3.5, loc = 1)
            axins.plot(dataset.time,  dataset.input_flows[flow_key])
            axins.scatter(dataset.time[t], dataset.input_flows[flow_key][t])
            axins.set_xlabel("time/ s")
            axins.set_ylabel("flow rate/ uL/h")

            plot_name = '{} time {}_line.png'.format(out_name,dataset.time[t])
            plt.savefig(plot_name, bbox_inches = "tight", transparent = True)
            image_names.append(plot_name)
            images.append(plt.imread(plot_name)) # read the written image into the images arra. imread is from imageio
            os.remove(plot_name)
            plt.close()
            print('frame') # just to see that it is working

    mimsave('{}.gif'.format(out_name), images, duration = 0.2) # imageio.mimsave converts the image array into a .gif. The duration kwarg is how long each frame lasts (in seconds, I think)
