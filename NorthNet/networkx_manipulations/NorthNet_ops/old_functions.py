'''
All functions below are old.
'''

def assign_rate_ranges(network, ranges):
    '''
    Needs updating
    '''
    for d in network.NetworkReactions:
        if network.NetworkReactions[d].Reaction_class == 'aldol_addition_no_stereo': # Finding aldol reactions and assigning them as the higher range if formaldehyde is in the reactants
            for r in network.NetworkReactions[d].Reactants:
                if r == 'C=O':
                    setattr(network.NetworkReactions[d], 'Range', ranges[network.NetworkReactions[d].Reaction_class]['High'])
                    break
                else:
                    setattr(network.NetworkReactions[d], 'Range', ranges[network.NetworkReactions[d].Reaction_class]['Low'])
                    break
        else:
            try:
                setattr(network.NetworkReactions[d], 'Range', ranges[network.NetworkReactions[d].Reaction_class]['High'])
            except:
                setattr(network.NetworkReactions[d], 'Range', ['unknown','unknown'])

def equations_from_network(network):
    '''
    Needs updating.
    '''
    '''Convert the network equations from SMILES reactions into a system of
    equations.'''

    species = [*network.NetworkCompounds]

    translation = {s:'' for s in species}
    j = 0
    for s in species:
        translation[s] = 'S[{}]'.format(j)
        j+=1

    reactions = [*network.NetworkReactions]

    rateconsts = {} # Create a list of rate constants (labels)
    ratemat = [['0' for x in species] for x in species]

    for l in range(0,len(reactions)):
        rateconsts[reactions[l]] = 'k[{}]'.format(l)

    for n in network.NetworkCompounds:

        ind1 = species.index(n)

        '''outgoing reactions'''
        for out in network.NetworkCompounds[n].Out:

            reacs = network.NetworkReactions[out].Reactants

            ki = rateconsts[out]
            if len(reacs) == 1:
                ratemat[ind1][ind1] +=  "-" + ki

            if len(reacs) == 2:
                z = reacs[:]
                z.remove(n)
                n2 = z[0]
                ind2 = species.index(n2)
                ratemat[ind1][ind2] +=  "-" + ki + '*' + translation[n]

            if len(reacs) == 3:
                z = reacs[:]
                z.remove(n)
                n2 = z[0]
                n3 = z[1]
                ind2 = species.index(n2)
                ratemat[ind1][ind2] +=  "-" + ki + '*' + translation[n] + '*' + translation[n3]

        for i in network.NetworkCompounds[n].In:
            reacs = network.NetworkReactions[i].Reactants
            ki = rateconsts[i]
            if len(reacs) == 1:
                n2 = reacs[0]
                ind2 = species.index(n2)
                ratemat[ind1][ind2] +=  '+' + ki

            if len(reacs) == 2:
                n2 = reacs[0]
                n3 = reacs[1]
                ind2 = species.index(n3)
                ratemat[ind1][ind2] +=  '+' + ki + '*' + translation[n2]

            if len(reacs) == 3:
                n2 = reacs[0]
                n3 = reacs[1]
                n4 = reacs[2]
                ind2 = species.index(n3)
                ratemat[ind1][ind2] +=  '+' + ki + '*' + translation[n2] + '*' + translation[n4]

    return ratemat,translation,rateconsts

def equations_from_network_classified_rates(network):
    '''
    Needs updating.
    '''
    '''Convert the network equations from SMILES reactions into a system of
    equations.'''

    species = [*network.NetworkCompounds]

    translation = {s:'' for s in species}
    j = 0
    for s in species:
        translation[s] = 'S[{}]'.format(j)
        j+=1

    reactions = [*network.NetworkReactions]

    rateconsts = {} # Create a list of rate constants (labels)
    ratemat = [['0' for x in species] for x in species]
    l = 0
    for e in network.NetworkReactions:
        if network.NetworkReactions[e].Reaction_class in [*rateconsts]:
            pass
        else:
            rateconsts[network.NetworkReactions[e].Reaction_class] = 'k[{}]'.format(l)
            l+=1

    for n in network.NetworkCompounds:
        ind1 = species.index(n)

        '''outgoing reactions'''
        for out in network.NetworkCompounds[n]['Out']:
            reacs = network.NetworkReactions[out].Reactants
            ki = rateconsts[network.NetworkReactions[out].Reaction_class]

            if len(reacs) == 1:
                ratemat[ind1][ind1] +=  "-" + ki

            elif len(reacs) == 2:
                z = reacs[:]
                z.remove(n)
                n2 = z[0]
                ind2 = species.index(n2)
                ratemat[ind1][ind2] +=  "-" + ki + '*' + translation[n]

            elif len(reacs) == 3:
                z = reacs[:]
                z.remove(n)
                n2 = z[0]
                n3 = z[1]
                ind2 = species.index(n2)
                ratemat[ind1][ind2] +=  "-" + ki + '*' + translation[n] + '*' + translation[n3]

        for i in network.NetworkCompounds[n]['In']:
            reacs = network.NetworkReactions[i].Reactants
            ki = rateconsts[network.NetworkReactions[i].Reaction_class]
            if len(reacs) == 1:
                n2 = reacs[0]
                ind2 = species.index(n2)
                ratemat[ind1][ind2] +=  '+' + ki

            if len(reacs) == 2:
                n2 = reacs[0]
                n3 = reacs[1]
                ind2 = species.index(n3)
                ratemat[ind1][ind2] +=  '+' + ki + '*' + translation[n2]

            if len(reacs) == 3:
                n2 = reacs[0]
                n3 = reacs[1]
                n4 = reacs[2]
                ind2 = species.index(n3)
                ratemat[ind1][ind2] +=  '+' + ki + '*' + translation[n2] + '*' + translation[n4]

    return ratemat,translation,rateconsts
