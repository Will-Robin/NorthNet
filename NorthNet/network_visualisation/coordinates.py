'''
Useful functions for manipulating networkx DiGraphs
'''

def get_network_lineplot(G):
    '''
    Create a numpy array from a DiGraph which can be used to plot its skeleton.

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
    Create a numpy array from a DiGraph which can be used to plot its skeleton
    as a scatter plot.

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
    Normalise the width and height of the DiGraph's coordinates.

    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    G2: networkx DiGraph
    '''
    import copy
    import numpy as np
    from NorthNet.network_visualisation import coordinates

    coords = coordinates.get_network_lineplot(G)

    net_width = (np.nanmax(coords[0])-np.nanmin(coords[0]))
    net_height = (np.nanmax(coords[1])-np.nanmin(coords[1]))

    G2 = copy.deepcopy(G)

    for n in G.nodes:
        pos = G.nodes[n]['pos']
        a = pos[0]/net_width
        b = pos[1]/net_height
        G2.nodes[n]['pos'] = (a,b)

    return G2

def rotate_network(G, radians):
    '''
    Rotate the DiGraph's coordinates

    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.
    radians: float
        rotation angle

    Returns
    -------
    G2: networkx DiGraph
    '''

    import numpy as np

    xy = get_network_scatter(G)

    ox, oy = xy[0].mean(),  xy[1].mean()

    G2 = copy.deepcopy(G)

    for n in G.nodes:
        px,py =  G.nodes[n]['pos']

        qx = ox + np.cos(radians) * (px - ox) - np.sin(radians) * (py - oy)
        qy = oy + np.sin(radians) * (px - ox) + np.cos(radians) * (py - oy)

        G2.nodes[n]['pos'] = (qx, qy)

    return G2
