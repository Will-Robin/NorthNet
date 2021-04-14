def score_reaction_types(reaction_network, scheme, amps_data, root_node):
    from NorthNet import network_generation as n_gen
    from NorthNet import conversions as conv

    score_def = lambda source,target: (source - target)/source

    type_scores = {reaction_network.NetworkReactions[r].ReactionTemplate.Name:0.0 for r in reaction_network.NetworkReactions}

    r_list = [reaction_network.NetworkReactions[r].ReactionSMILES for r in scheme]

    r_net = n_gen.network_from_reaction_list(r_list,network_name = "")
    G = conv.convert_to_networkx(r_net)

    ord_amps =  [x for v,x in sorted( zip( list(amps_data.values()), [*amps_data] ) ) if x != "C=O"]
    ord_amps.reverse()

    if root_node in ord_amps:
        ord_amps.remove(root_node)
        ord_amps.insert(0,root_node)

    for x in range( 0, len(ord_amps) - 1 ):
        for y in range( x+1, len(ord_amps)):
            loc_rxns = []

            amp_drop = score_def(amps_data[ord_amps[x]], amps_data[ord_amps[y]])

            if nx.has_path(G, ord_amps[x], ord_amps[y]):
                path = nx.shortest_path(G, source = ord_amps[x], target = ord_amps[y])
                for p in path:
                    if ">>" in p:
                        loc_rxns.append(reaction_network.NetworkReactions[p].ReactionTemplate.Name)
                loc_rxns = list(set(loc_rxns))
                for l in loc_rxns:
                    type_scores[l] += amp_drop/len(loc_rxns) # bigget relative amplitude drop = higher degreee of rate control?

    return type_scores

def compare_network_pathways(net1,net2):
    '''
    net1: NorthNet Network object
        Parent network
    net2: NorthNet Network object
        Child network
    '''

    # 1. try to compare relative number of reaction classes
    parent_class_counts = {}
    for r in net1.NetworkReactions:
        r_class = net1.NetworkReactions[r].ReactionTemplate.Name
        if  r_class in parent_class_counts:
            parent_class_counts[r_class] += 1
        else:
            parent_class_counts[r_class] = 1

    child_class_counts = {k:0 for k in parent_class_counts}
    for r in net2.NetworkReactions:
        r_class = net2.NetworkReactions[r].ReactionTemplate.Name
        child_class_counts[r_class] += 1

    comparisons = {}
    for p in parent_class_counts:
        comparisons[p] = child_class_counts[p]/parent_class_counts[p]


    comparisons = {k:comparisons[k] for k in sorted([*comparisons], key = lambda x:len(x))}
    return comparisons
