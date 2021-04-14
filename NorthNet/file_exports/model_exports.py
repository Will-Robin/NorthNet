def write_model_as_module(network, numba_decoration = False):

    get_index = lambda x: int(x[x.find("[")+1:x.find("]")])
    eq_text = network_ops.write_model_equation_text(network).split("\n")
    #mat_text = network_ops.write_model_matrix_text(network)
    jac_text = network_ops.write_Jacobian_matrix_text(network)

    species_indices, rate_consts_indices, inflows = network_ops.network_indices(network)

    fname = "{}.py".format(network.Name)

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

def write_model_to_file(filename,ratemat,network,translation,rateconsts, generation_protocol):
    '''
    Needs updating.
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
