def set_network_coords(network, pos):
    import numpy as np

    xmin = np.mean([pos[p][0] for p in pos])
    ymin = np.mean([pos[p][1] for p in pos])
    for n in network.nodes:
        if n in pos:
            network.nodes[n]['pos'] = pos[n]
        else:
            network.nodes[n]['pos'] = (xmin, ymin)

    return network

def get_network_lineplot(G):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    net_lines: numpy 2D array
        Coordinates for plotting a line plot of the network.
    '''

    import numpy as np

    net_lines = []
    for e in G.edges:
        for n in e:
            net_lines.append(G.nodes[n]["pos"])
        net_lines.append((np.nan,np.nan))

    net_lines = np.array(net_lines)
    net_lines = net_lines.T
    return net_lines

def get_network_scatter(G):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    net_lines: numpy 2D array
        Coordinates for plotting a line plot of the network.
    '''

    import numpy as np

    net_scatter = []
    for n in G.nodes:
        net_scatter.append(G.nodes[n]["pos"])

    net_scatter = np.array(net_scatter)
    net_scatter = net_scatter.T
    return net_scatter


def normalise_network_coordinates(G):
    '''
    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    None
    '''
    import numpy as np
    from NorthNet.network_manipulations.networkx_ops import coordinates
    coords = coordinates.get_network_lineplot(G)

    net_width = (np.nanmax(coords[0])-np.nanmin(coords[0]))
    net_height = (np.nanmax(coords[1])-np.nanmin(coords[1]))

    for n in G.nodes:
        pos = G.nodes[n]['pos']
        a = pos[0]/net_width
        b = pos[1]/net_height
        G.nodes[n]['pos'] = (a,b)

def rotate_network(G, radians, origin):
    import numpy as np

    xy = get_network_scatter(G)
    x, y = xy[0], xy[1]

    offset_x, offset_y = origin
    adjusted_x = (x - offset_x)
    adjusted_y = (y - offset_y)
    cos_rad = np.cos(radians)
    sin_rad = np.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y

    for c,n in enumerate(G.nodes):
        G.nodes[n]['pos'] = (qx[c], qy[c])


def set_network_coordinates(G, coords_file):
    '''
    Adds coordinates to node attributes from a .csv file.

    Parameters
    ----------
    G: networkx DiGraph
        Network
    coords_file: str
        Path to file containing node coordinates (format: compound, x, y newline)

    '''
    # Build coordinates list
    spec_coords = {}
    with open(coords_file, "r") as f:
        for line in f:
            ln = line.strip("\n")
            ln = ln.split(",")
            if ln[0].strip('"') in spec_coords:
                pass
            else:
                spec_coords[ln[0].strip('"')] = tuple([float(x) for x in ln[1:]])

    for n in G.nodes:
        if n in spec_coords:
            G.nodes[n]["pos"] = spec_coords[n]
        else:
            G.nodes[n]["pos"] = (randint(0,100),randint(0,100))

    return G
