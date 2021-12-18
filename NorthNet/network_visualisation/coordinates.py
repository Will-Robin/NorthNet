import copy
import numpy as np
from NorthNet.network_visualisation import coordinates
'''
Useful functions for manipulating networkx DiGraphs
'''
def set_network_coords(G, pos):
    '''
    Add coordinate information into a networkx DiGraph.

    Parameters
    ----------
    G: networkx DiGraph
        Graph to extract nodes from.

    Returns
    -------
    net_lines: numpy 2D array
        Coordinates for plotting a line plot of the network.
    '''


    xmin = np.mean([pos[p][0] for p in pos])
    ymin = np.mean([pos[p][1] for p in pos])
    for node in G.nodes:
        if node in pos:
            G.nodes[node]['pos'] = pos[node]
        else:
            G.nodes[node]['pos'] = (xmin, ymin)

    return G

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


    
    net_lines = []
    for edge in G.edges:
        for node in edge:
            net_lines.append(G.nodes[node]["pos"])
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

    net_scatter = []
    for node in G.nodes:
        net_scatter.append(G.nodes[node]["pos"])

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

    coords = coordinates.get_network_lineplot(G)

    net_width = (np.nanmax(coords[0])-np.nanmin(coords[0]))
    net_height = (np.nanmax(coords[1])-np.nanmin(coords[1]))

    G2 = copy.deepcopy(G)

    for node in G.nodes:
        pos = G.nodes[node]['pos']
        x_coordinate = pos[0]/net_width
        y_coordinate = pos[1]/net_height
        G2.nodes[node]['pos'] = (x_coordinate,y_coordinate)

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

    xy = get_network_scatter(G)

    ox, oy = xy[0].mean(),  xy[1].mean()

    G2 = copy.deepcopy(G)

    for node in G.nodes:
        px,py =  G.nodes[node]['pos']

        qx = ox + np.cos(radians) * (px - ox) - np.sin(radians) * (py - oy)
        qy = oy + np.sin(radians) * (px - ox) + np.cos(radians) * (py - oy)

        G2.nodes[node]['pos'] = (qx, qy)

    return G2

