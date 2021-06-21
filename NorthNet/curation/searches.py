def get_hits_single_criterion(reaction_list, attribute, criterion):
    '''
    Parameters
    ----------
    Returns
    -------
    '''

    tag_list = [] # list which will contain boolean tags which will allow items to be eliminated based on a condition
    for r in reaction_list:
        tag = False
        if criterion in r.Info[attribute]: # if the element is in the reaction string, the reaction is tagged for removal
            tag = True
        tag_list.append(tag)

    test_fails = [x for x,y in zip(reaction_list,tag_list) if y == False]
    for t in test_fails:
        r.Info["Rejection_reason"] = "{} did not satify {}.".format(attribute, criterion)

    reaction_list = [x for x,y in zip(reaction_list,tag_list) if y == True]

    return reaction_list, test_fails
