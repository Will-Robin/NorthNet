def remove_on_list_criterion(reaction_list, attribute, out_list):
    '''
    Parameters
    ----------
    reaction_list: list
        List of reaction objects.
    attribute: str
        Attribute on which selection will be based.
    out_list: list
        List of
    Returns
    -------
    '''

    tag_list = [] # list which will contain boolean tags which will allow items to be eliminated based on a condition
    for r in reaction_list:
        tag = False
        for o in out_list: # cycle through the precious/toxic metals list
            if o in r.Info[attribute]: # if the element is in the reaction string, the reaction is tagged for removal
                tag = True
                r.Info["Rejection_reason"] = "Rejected: contains {} (field: {})".format(o, attribute)
        tag_list.append(tag)

    test_fails = [x for x,y in zip(reaction_list,tag_list) if y == True]
    reaction_list = [x for x,y in zip(reaction_list,tag_list) if y == False]

    return reaction_list, test_fails

def remove_empty_reaction_SMILES(reaction_list):
    '''
    Parameters
    ----------
    Returns
    -------
    '''

    tag_list = []
    for rx in reaction_list:
        tag = False
        if len(rx.ReactionSMILES) == 0: # find empty reaction strings and remove them.
            tag = True
        tag_list.append(tag)

    test_fails = [x for x,y in zip(reaction_list,tag_list) if y == True]

    for t in test_fails:
        t.Info["Rejection_reason"] = "Rejected: Empty SMILES"

    reaction_list = [x for x,y in zip(reaction_list,tag_list) if y == False]

    return reaction_list, test_fails

def reaction_string_replacements(reaction_list, replacements_dict):
    '''
    Parameters
    ----------
    Returns
    -------
    '''

    for rx in reaction_list:
        print("Performing replacements in", rx)

        for k in replacements_dict.keys():
            rx.ReactionSMILES =  rx.ReactionSMILES.replace(k,replacements_dict[k])

    return reaction_list

def remove_individual_compounds(reaction_list, remove_list):
    '''
    Parameters
    ----------
    Returns
    -------
    '''

    for rx in reaction_list:
        reactants = rx.ReactionSMILES.split(">>")[0].split(".")
        products = rx.ReactionSMILES.split(">>")[1].split(".")
        reactants = [r for r in reactants if r not in remove_list]
        products  = [p for p in products if p not in remove_list]

        rx.ReactionSMILES = reaction_ops.reconstitute_reaction(reactants,products)

    return reaction_list

def metal_ion_dissociation(reaction_list, metals_list):
    '''
    Parameters
    ----------
    Returns
    -------
    '''
    print("metal_ion_dissociation function needs fixing before use. e.g. update reaction and new mol attributes")
    for r in reaction_list:
        for m in metals_list:
            if m in r.ReactionSMILES:
                replacement = r.ReactionSMILES.replace(m,'.')

                reac_string = ''
                prod_string = ''

                reactants = rx.ReactionSMILES.split(">>")[0].split(".")
                products = rx.ReactionSMILES.split(">>")[1].split(".")
                tag = False

                for reac in reactants:
                    if Chem.MolFromSmiles(reac) == None:
                        pass
                    else:
                        reac_string+=reac + '.'

                reac_string = reac_string[:-1]

                for prod in products:
                    if Chem.MolFromSmiles(prod) == None:
                        pass
                    else:
                        prod_string+=reac + '.'

                prod_string = prod_string[:-1]

                replacement = reac_string + ">>" + prod_string

    return reaction_list

def canonicalise_reactions(reaction_list):
    '''

    This should be taken care of upon initisation of the reaction objects,
    but performed again to take into account string manipulations during curation.

    Parameters
    ----------
    Returns
    -------
    '''

    for rx in reaction_list:
        reactants = rx.ReactionSMILES.split(">>")[0].split(".")
        products = rx.ReactionSMILES.split(">>")[1].split(".")

        reactants = [molecule_ops.canonicalise(r) for r in reactants]
        products  = [molecule_ops.canonicalise(p) for p in products]

        rx.ReactionSMILES = reaction_ops.reconstitute_reaction(reactants,products)

    return reaction_list

def remove_self_self_reactions(reaction_list):
    '''
    Parameters
    ----------
    Returns
    -------
    '''

    tag_list = []
    for rx in reaction_list:
        reactants = rx.ReactionSMILES.split(">>")[0].split(".")
        products = rx.ReactionSMILES.split(">>")[1].split(".")
        tag = False
        for r in reactants:
            if any(item in r for item in ['C', 'c']):
                for p in products:
                    if r == p:
                        tag = True
        tag_list.append(tag)

    test_fails = [x for x,y in zip(reaction_list,tag_list) if y == True]
    for t in test_fails:
        t.Info["Rejection_reason"] = "Rejected: Self>>Self"

    reaction_list = [x for x,y in zip(reaction_list,tag_list) if y == False]

    return reaction_list, test_fails

def remove_repeated_components(reaction_list):
    '''
    Parameters
    ----------
    Returns
    -------
    '''


    for r in reaction_list:
        reactants = reaction_object.ReactionSMILES.split(">>")[0].split(".")
        products = reaction_object.ReactionSMILES.split(">>")[1].split(".")

        for x in range(0,len(reactants)-1):

            mols = [Chem.MolFromSmiles(reactants[x]),Chem.MolFromSmiles(reactants[x+1])]
            fps  = [Chem.RDKFingerprint(m) for m in mols]
            sim = DataStructs.FingerprintSimilarity(fps[0],fps[1])

            if sim == 1.0:
                reactants = reactants[x:]
            else:
                pass

        for x in range(0,len(r.ProductsAsMols)-1):
            mols = [r.ProductsAsMols[x],r.ProductsAsMols[x+1]]
            fps  = [Chem.RDKFingerprint(m) for m in mols]
            sim = DataStructs.FingerprintSimilarity(fps[0],fps[1])

            if sim == 1.0:
                products = products[x:]
            else:
                pass

        r.ReactionSMILES = reaction_ops.reconstitute_reaction(reactants,products)
