import networkx as nx
from scipy.cluster import hierarchy

def graph_from_linkage(linkage_mat, id_modifier = ''):
    '''
    linkage_mat: scipy condensed linkage matrix
        linkage matrix

    id_modifier: str
        name modifier (can be used to distinguish different clusterings)
    '''
    rootnode, nodelist = hierarchy.to_tree(linkage_mat, rd= True)

    G = nx.DiGraph()

    for n in nodelist:
        source = '{}{}'.format(n.id, id_modifier)
        G.add_node(source, distance = n.dist, leaf = n.is_leaf(),
                    id = n.id)
        if n.count > 1:
            target_left = '{}{}'.format(n.left.id, id_modifier)

            G.add_node(target_left, distance = n.left.dist, leaf = n.left.is_leaf(),
                        id = n.left.id)
            G.add_edge(source, target_left, weight = abs(n.left.dist - n.dist),
                        dist = abs(n.left.dist - n.dist))

            target_right = '{}{}'.format(n.right.id, id_modifier)
            G.add_node(target_right, distance = n.right.dist, leaf = n.right.is_leaf(),
                        id = n.right.id)
            G.add_edge(source, target_right, weight = abs(n.right.dist - n.dist),
                        dist = abs(n.right.dist - n.dist))

    return G
