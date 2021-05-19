def map_network_reactions(network, env_radius = 1):
    '''
    Peform atom-atom mapping for all reactions in a network.
    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network containing reactions to be mapped.
    Returns
    -------
    None
    '''
    # Map reactions
    for r in network.NetworkReactions:
        mapped_r = r_o.aam_smiles(network.NetworkReactions[r].ReactionSMILES)
        sub_reaction = r_o.extract_reacting_substructures(mapped_r, radius = env_radius)

        network.NetworkReactions[r].ReactionTemplate = Classes.Reaction_Template(sub_reaction,sub_reaction,sub_reaction.split(">>")[0].split("."), sub_reaction.split(">>")[1].split("."))
        network.NetworkReactions[r].MappedReaction = mapped_r

        for d in network.NetworkReactions[r].Database_Entries:
            d.Info["mapped_reaction"] = mapped_r
            d.Info["sub_reaction"] = sub_reaction

def extract_reaction_centres(network, env_radius = 1):
    '''
    Peform atom-atom mapping for all reactions in a network.
    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network containing reactions to be mapped.
        Requires mapped reactions.
    Returns
    -------
    None
    '''
    pass
