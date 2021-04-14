def scheme_dict_to_network_dict(schemes):
    from NorthNet import network_generation as n_gen
    from NorthNet import conversions as conv
    rxn_networks = {}
    for s in schemes:
        net = n_gen.network_from_reaction_list(schemes[s],network_name = "{}".format(s))
        G = conv.convert_to_networkx(net)
        for n in ["C=O", 'O', '[OH-]']:
            if n in G.nodes:
                G.remove_node(n)
        rxn_networks[s] = G
    return rxn_networks
