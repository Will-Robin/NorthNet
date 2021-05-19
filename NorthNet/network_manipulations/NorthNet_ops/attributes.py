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
