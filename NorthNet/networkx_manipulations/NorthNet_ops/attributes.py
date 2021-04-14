def network_indices(network):

    compounds = [x for x in network.NetworkCompounds]
    reactions = [*network.NetworkReactions]

    species = {s:"S[{}]".format(c) for c,s in enumerate(compounds)}
    rate_consts = {k:"k[{}]".format(c) for c,k in enumerate(reactions)}

    out_ind = -1
    for c,k in enumerate(rate_consts):
        if k.endswith(">>") and out_ind == -1:
            out_ind = c
        if k.endswith(">>") and out_ind != -1:
            rate_consts[k] = "k[{}]".format(out_ind)

    inflows = {k:"C[{}]".format(c) for c,k in enumerate([x for x in network.NetworkReactions if x.startswith(">>")])}

    return species, rate_consts, inflows

def sort_mols(network):
    '''
    Sorts molecules in a network by carbon.

    Will be replaced.

    Parameters
    ----------
    network: NorthNet network object
        Network containg molecules to be sorted.
    Returns
    -------
    sort: list
        A sorted list of SMILES.
    '''
    print("Consider using method in the Reaction class instead.")
    '''Sorts molecules in network list by carbon number.'''
    ind = [] # empty list for carbon numbers by which the list will be sorted
    for m in network:
        cut = m.count("C")+m.count('c') # get number of carbons in species
        ind.append(cut) # add the number of carbons to a list (effectively indexed to the network)
    sort = [x for _,x in sorted(zip(ind,network))] # sort the network based on the ind list
    return sort
