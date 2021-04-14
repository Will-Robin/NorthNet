def write_references(network,fname):
    '''
    Places the reactions in a network into a spreadsheet (comma separated values)
    next to their literature references.
    Parameters
    ----------
    network: NorthNet Network object
        Network to be converted to spreadheet.
    fname: str
        A name for the output file.
    Returns
    -------
    None
    '''
    with open("{}_references.csv".format(fname), "w") as f:
        for r in network.NetworkReactions:
            f.write("{},".format(r))
            for d in network.NetworkReactions[r].Database_Entries:
                f.write("{},".format(d.Info["References"].replace(","," ")))
            f.write("\n")

def write_assigned_equation_list(filename, network, rateconsts, generation_protocol):
    '''
    Needs updating

    Parameters
    ----------
    filename: str
        pass
    network: NorthNet Network Object
        pass
    rateconsts: ???
        pass
    generation_protocol: ???
        pass
    Returns
    -------
    None
    '''
    with open('{}_chemical_equation_list.csv'.format(filename), 'w') as f: # writing the data to a .csv file.

        f.write('generation_protocol:\n')
        for g in generation_protocol:
            reagents = '+'.join(g[1])
            f.write("{} + {},".format(g[0],reagents))
        f.write("\n")

        f.write('reaction,type,lower rate limit/ s-1 or M-1s-1,upper rate limit/ s-1 or M-1s-1,assignment')
        f.write('\n')
        for r in [*network.NetworkReactions]:
            ins  = int( rateconsts[network.NetworkReactions[r].Generation_Details[0]].replace('[','').replace(']','').replace('k',''))
            f.write('{},{},{},{},{},'.format(r,network.NetworkReactions[r].Generation_Details[0],network.NetworkReactions[r][0].Range[0],network.NetworkReactions[r][0].Range[1],ins))
            f.write('\n')

def write_equation_list(filename, network, rateconsts, generation_protocol):
    '''
    Needs updating
    '''
    with open('{}_chemical_equation_list.csv'.format(filename), 'w') as f: # writing the data to a .csv file.
        f.write('generation_protocol:\n')
        for g in generation_protocol:
            reagents = '+'.join(g[1])
            f.write("{} + {},".format(g[0],reagents))
        f.write("\n")

        f.write('reaction,type,lower rate limit/ s-1 or M-1s-1,upper rate limit/ s-1 or M-1s-1,assignment')
        f.write('\n')
        for r in [*network.NetworkReactions]:
            ins  = int(rateconsts[r].replace('[','').replace(']','').replace('k',''))
            f.write('{},{},{},{},{},'.format(r,network.NetworkReactions[r].Reaction_class,network.NetworkReactions[r].Range[0],network.NetworkReactions[r].Range[1],ins))
            f.write('\n')

    return 'This is for writing a file!'
