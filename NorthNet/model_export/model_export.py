from NorthNet.model_export import model_export

def network_indices(network):
    '''
    Parameters
    ----------
    network: NorthNet ReactionNetwork object
        Network for which model will be written.

    Returns
    -------
    species, rate_consts, inflows
        Dictionaries of tokens and their indices.
    '''

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

    inflows = {}
    inflow_reactions = [x for x in network.NetworkReactions if x.startswith(">>")]
    for c,k in enumerate(inflow_reactions):
        inflows[k] = "C[{}]".format(c)

    return species, rate_consts, inflows

def write_model_as_module(network, numba_decoration = False, filename = ''):
    '''
    Parameters
    ----------
    network: NorthNet ReactionNetwork object
        Network for which model will be written.

    numba_decoration: bool
        Whether to add numba decorator (wrapper function).
    filename: str
        Name for file.
    '''

    get_index = lambda x: int(x[x.find("[")+1:x.find("]")])

    eq_text = model_export.write_model_equation_text(network).split("\n")
    # mat_text = model_export.write_model_matrix_text(network)
    jac_text = model_export.write_Jacobian_matrix_text(network)

    species_indices, rate_consts_indices, inflows = model_export.network_indices(network)

    if filename == '':
        fname = "{}.py".format(network.Name)
    else:
        fname = filename

    with open(fname, "w") as f:
        f.write("import numpy as np\n")
        if numba_decoration:
            f.write("import numba\n")
        f.write("\n")
        if numba_decoration:
            f.write("@numba.jit(numba.float64[:](numba.float64,numba.float64[:],numba.float64[:],numba.float64[:]),\n"
                       "\tlocals={'P': numba.float64[:],'k': numba.float64[:]},nopython=True)\n")
        f.write("def model_function(time, S, k_base, C):\n")
        f.write("\n")
        f.write("\tP = np.zeros(len(S))\n")
        f.write("\n")
        f.write("\tk = np.copy(k_base)\n")
        f.write("\n")
        for e in eq_text:
            f.write("\t{}\n".format(e))
        f.write("\treturn P\n")
        f.write("\n")

        f.write("def wrapper_function(time, S, k_base, C):\n")
        f.write("\treturn model_function(time, S, k_base, C)\n")
        f.write("\n")

        f.write("species = {")
        [f.write("'{}':{},".format(k,get_index(species_indices[k]))) for k in species_indices]
        f.write("}\n")
        f.write("\n")

        f.write("reactions = {")
        [f.write("'{}':{},".format(k,get_index(rate_consts_indices[k]))) for k in rate_consts_indices]
        f.write("}\n")
        f.write("\n")

        f.write("inputs = {")
        [f.write("'{}':{},".format(k,get_index(inflows[k]))) for k in inflows]
        f.write("}\n")
        f.write("\n")

        f.write("k = np.zeros(max(reactions.values())+1) # rate constants\n")
        f.write("\n")

        f.write("S = np.zeros(len(species)) # initial concentrations\n")
        f.write("\n")

        f.write("C = np.zeros(len(inputs)) # input concentrations\n")
        f.write("\n")

def write_model_equation_text(network):
    '''
    Parameters
    ----------
    network: NorthNet ReactionNetwork object
        Network to be written.

    Returns
    -------
    eq_text: str
        Rate equations in text form.
    '''

    compounds = [x for x in network.NetworkCompounds]
    reactions = [*network.NetworkReactions]

    species, rate_consts, inflows = model_export.network_indices(network)

    eq_text = ""

    for count,c in enumerate(compounds):
        eq_text += "P[{}] = ".format(count)
        for i in network.NetworkCompounds[c].In:
            ki = "+{}*".format(rate_consts[i])

            if len(network.NetworkReactions[i].Reactants) == 0:
                specs = inflows[i]
            else:
                specs = "*".join([species[x] for x in network.NetworkReactions[i].Reactants])

            eq_text += "{}{}".format(ki,specs)

        for o in network.NetworkCompounds[c].Out:
            ki = "-{}*".format(rate_consts[o])
            specs = "*".join([species[x] for x in network.NetworkReactions[o].Reactants])
            eq_text += "{}{}".format(ki,specs)

        eq_text += "\n"

    return eq_text

def write_model_matrix_text(network):
    '''
    Parameters
    ----------
    network: NorthNet ReactionNetwork object
        Network to be written.

    Returns
    -------
    mat_text: str
        Rate equations in numpy matrix form.
    '''
    compounds = [x for x in network.NetworkCompounds]
    reactions = [*network.NetworkReactions]

    species, rate_consts, inflows = model_export.network_indices(network)

    ratemat = [['0' for x in species] for x in species]

    for c in compounds:

        ind1 = compounds.index(c)

        '''outgoing reactions'''
        for i in network.NetworkCompounds[c].In:
            reacs = network.NetworkReactions[i].Reactants
            ki = rate_consts[i]

            if len(reacs) == 0:
                ratemat[ind1][ind1] +=  "(+{}*{})/{}".format(ki,inflows[c],species[c])

            if len(reacs) == 1:
                n2 = reacs[0]
                ind2 = compounds.index(n2)
                ratemat[ind1][ind2] +=  '+' + ki

            if len(reacs) == 2:
                n2 = reacs[0]
                n3 = reacs[1]
                ind2 = compounds.index(n3)
                ratemat[ind1][ind2] +=  '+' + ki + '*' + species[n2]

            if len(reacs) == 3:
                n2 = reacs[0]
                n3 = reacs[1]
                n4 = reacs[2]
                ind2 = compounds.index(n3)
                ratemat[ind1][ind2] +=  '+' + ki + '*' + species[n2] + '*' + species[n4]

        for out in network.NetworkCompounds[c].Out:

            reacs = network.NetworkReactions[out].Reactants

            ki = rate_consts[out]
            if len(reacs) == 1:
                ratemat[ind1][ind1] +=  "-" + ki

            if len(reacs) == 2:
                z = reacs[:]
                z.remove(c)
                n2 = z[0]
                ind2 = compounds.index(n2)
                ratemat[ind1][ind2] +=  "-" + ki + '*' + species[c]

            if len(reacs) == 3:
                z = reacs[:]
                z.remove(c)
                n2 = z[0]
                n3 = z[1]
                ind2 = compounds.index(n2)
                ratemat[ind1][ind2] +=  "-" + ki + '*' + species[c] + '*' + species[n3]

    mat_text = "["
    for r in ratemat:
        mat_text += "[" + ",".join(r) + "],\n"

    mat_text = mat_text.strip(",\n") + "]"

    return mat_text

def write_Jacobian_matrix_text(network):
    '''
    Parameters
    ----------
    network: NorthNet ReactionNetwork object
        Network to be written.

    Returns
    -------
    jac_text: str
        Jacobian matrix as text.
    '''
    compounds = [x for x in network.NetworkCompounds]
    reactions = [*network.NetworkReactions]

    species, rate_consts, inflows = model_export.network_indices(network)

    jac_mat = [['0' for x in species] for x in species]

    for c,comp1 in enumerate(compounds):
        for c2,comp2 in enumerate(compounds):
            element = ""
            for i in network.NetworkCompounds[comp1].In:
                if comp2 in network.NetworkReactions[i].Reactants:
                    reacs = [species[x] for x in network.NetworkReactions[i].Reactants if x != comp2]
                    ki = "+{}".format(rate_consts[i])
                    element += "{}*{}".format(ki,"*".join(reacs))
                else:
                    pass

            for o in network.NetworkCompounds[comp1].Out:
                if comp2 in network.NetworkReactions[o].Reactants:
                    reacs = [species[x] for x in network.NetworkReactions[o].Reactants if x != comp2]
                    ki = "-{}".format(rate_consts[o])
                    element += "{}*{}".format(ki,"*".join(reacs))
                else:
                    pass

            jac_mat[c][c2] += element

    jac_text = ""
    for r in jac_mat:
        jac_text += "[" + ",".join(r) + "],\n"

    jac_text = jac_text.strip(",\n") + "]"

    return jac_text

def write_model_to_file(filename,ratemat,network,translation,rateconsts, generation_protocol):
    '''
    Not such a good framework. Needs updating.
    '''
    '''
    Writes model output
    Parameters
    ----------
    filename: str
        pass
    ratemat: ??
        pass
    network: NorthNet Network Object
        pass
    translation: ???
        pass
    rateconsts: ???
        pass
    generation_protocol:??
        pass

    Returns
    -------
    None
    '''
    with open('{}_model_topology.txt'.format(filename), 'w') as f:
        f.write("Model"+"\n")
        f.write("{}\n".format(network.Name))

        f.write("rate_matrix_start\n")
        for x in ratemat:
            l = map(str,x)
            f.write("," + "\n")
            f.write("["+",".join(l) + "]")
        f.write('\n')
        f.write("rate_matrix_end\n")
        f.write('\n')

        v = map(str,translation.values())

        f.write("species_symbols_start\n")
        f.write(",".join(v) + "\n")
        f.write("species_symbols_end:\n")
        f.write('\n')

        f.write("species_identities_start\n")
        f.write('{')
        for tr in translation:
            f.write("'{}':{},".format(tr, translation[tr].replace("[","").replace("]","").replace("S","")))
        f.write('}\n')
        f.write("species_identities_end\n")
        f.write('\n')

        k_array = map(str,rateconsts.values())
        f.write("rate_constant_symbols_start\n")
        f.write(",".join(k_array) + "\n")
        f.write("rate_constant_symbols_end\n")
        f.write('\n')

        f.write("rate_constant_assignments_start\n")
        f.write('{')
        for tr in rateconsts:
            f.write("'{}':{},".format(tr,rateconsts[tr].replace("[","").replace("]","").replace("k","")))
        f.write('}\n')
        f.write("rate_constant_assignments_end\n")
        f.write('\n')

        f.write('generation_protocol:\n')
        for g in generation_protocol:
            reagents = '+'.join(g[1])
            f.write("{} + {},".format(g[0],reagents))
        f.write("\n")
