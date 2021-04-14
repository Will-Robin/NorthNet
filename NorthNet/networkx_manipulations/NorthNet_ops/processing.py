def reaction_temp_process(network):
    '''
    Processes the Temperature attributes the reactions in a network.
    If two values or a range are given, the temperatuer is averaged.
    If there is no temperature given, the temp is set to 25 oC.
    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network with temperature information inside it.
    Returns
    -------
    None
    '''
    temp_key = "Temperature (Reaction Details) [C]"
    for r in network.NetworkReactions:
        for x in network.NetworkReactions[r].Database_Entries:
            if ";" in x.Info[temp_key]:
                t = x.Info[temp_key].split(";")
                if "- " in t[0]:
                    spl = t[0].split("- ")
                    x.Info[temp_key] = (float(spl[0]) + float(spl[1]))/2
                else:
                    x.Info[temp_key] = float(t[0])
            elif " - " in x.Info[temp_key]:
                t = x.Info[temp_key].split(" - ")
                x.Info[temp_key] = (float(t[0])+float(t[1]))/2
            elif x.Info[temp_key] == "":
                x.Info[temp_key] = 25.0
            else:
                x.Info[temp_key] = float(x.Info[temp_key])

def balance_reaction_network(network):
    '''
    Balancing the stiochiometry of all reactions in a network.
    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be balanced.

    Returns
    -------
    None
    '''
    # some quick cleaning
    elems = ["C","N","S","O","Br","Cl","I","F","P"]
    for r in network.NetworkReactions:
        for e in elems:
            r_o.balance_reaction(network.NetworkReactions[r], e)
