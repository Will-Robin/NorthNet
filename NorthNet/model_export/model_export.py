from NorthNet.model_export import model_export
from NorthNet import info_params
from NorthNet import Classes
import numpy as np

def get_flow_rate(flow_profiles, time_limit = 1e100):

    '''
    Parameters
    ----------
    flow_profiles: dict
        dictionary of flow profiles

    time_limit: float
        Max time to include in the flow profile.

    Returns
    -------
    total_flows: 1D numpy array
        Numpy arrays total flow over time

    '''

    conc_flow_key_pairs = []
    for fl in flow_profiles:
        if 'time' in fl:
            time_axis = flow_profiles[fl]
        elif '/ M' in fl:
            inp_name = fl.strip('/ M')
            for fl2 in [*flow_profiles]:
                if inp_name in fl2 and fl2 != fl:
                    conc_flow_key_pairs.append((fl,fl2))

    idx = np.where(time_axis < time_limit)[0]
    total_flows = np.zeros(len(time_axis))
    for fl in flow_profiles:
        if 'flow' in fl and not 'time' in fl:
            total_flows += flow_profiles[fl]

    # ASSUMING THAT THE FLOW RATE IS IN UNITS OF uL/h
    # coverting to L/s
    total_flows /= 1e6
    total_flows /= 3600
    # convert to residence time
    total_flows = total_flows/(411*1e-6)

    return time_axis[idx], total_flows[idx]

def concentrations_from_flow_profile(flow_profiles, time_limit = 1e100):

    '''
    Parameters
    ----------
    flow_profiles: dict
        dictionary of flow profiles

    time_limit: float
        Max time to include in the flow profile.

    Returns
    -------
    concentrations: dict
        Numpy arrays of concentration inputs over time.
    '''

    conc_flow_key_pairs = []
    for fl in flow_profiles:
        if 'time' in fl:
            time_axis = flow_profiles[fl]
        elif '/ M' in fl:
            inp_name = fl.strip('/ M')
            for fl2 in [*flow_profiles]:
                if inp_name in fl2 and fl2 != fl:
                    conc_flow_key_pairs.append((fl,fl2))

    total_flows = np.zeros(len(time_axis))
    for fl in flow_profiles:
        if 'flow' in fl and not 'time' in fl:
            total_flows += flow_profiles[fl]

    idx = np.where(time_axis < time_limit)[0]
    concentrations = {}
    for p in conc_flow_key_pairs:
        moles = flow_profiles[p[0]]*flow_profiles[p[1]]#/total_flows
        conc = moles/(411*1e-6)
        if 'NaOH' in p[0]:
            concentrations['[OH-]/ M'] = conc[idx]
        else:
            concentrations[p[0]] = conc[idx]

    return time_axis[idx], concentrations

def write_flow_profile_text(flow_profiles, network, time_limit = 1e100):
    '''
    Parameters
    ----------
    flow_profiles: dict
        flow profiles
    network: NorthNet Network objects
        Network containing flow_terms

    time_limit: float
        Time limit for flow profile calculation.

    Returns
    -------
    text: numpy arrays for inputs

    '''
    get_index = lambda x: int(x[x.find("[")+1:x.find("]")])

    t, total_flow = model_export.get_flow_rate(flow_profiles, time_limit = time_limit)

    species, rate_consts, inflows, flow_ins, flow_outs = model_export.network_indices(network)


    collection_array = np.zeros((len(conc_dict)+2, len(t)))
    collection_array[0] = t
    for i in inflows:
        if 'O_#0' in i:
            pass
        else:
            key = i.strip('_#0') + '/ M'
            key = key.replace(">>","")
            idx = get_index(inflows[i])
            collection_array[idx+1] = conc_dict[key]

    collection_array[-1] = total_flow

    text = 'allF = np.array('
    text += np.array2string(collection_array,
                                 formatter={'float_kind':lambda x: "%.9f" % x},
                                 separator=',',threshold=np.inf)
    text += ')'

    return text


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

    species, rate_consts, inflows, flow_ins, flow_outs = model_export.network_indices(network)

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

def write_model_as_module_text_B(network, flow_profiles,
                                 time_limit = 1e100,
                                 numba_decoration = False):
    '''
    Parameters
    ----------
    network: NorthNet Network object
        Network for which model will be written.
    flow_profiles: dict
        flow profiles
    time_limit: float
        Time limit for flow profile calculation.
    numba_decoration: bool
        Whether to add numba decorator (wrapper function).
    filename: str
        Name for file.
    '''

    get_index = lambda x: int(x[x.find("[")+1:x.find("]")])

    eq_text = model_export.write_model_equation_text(network).split("\n")
    # mat_text = model_export.write_model_matrix_text(network)
    jac_text = model_export.write_Jacobian_matrix_text(network)

    flow_profile_text = model_export.write_flow_profile_text(flow_profiles, network,
                                                time_limit = time_limit)

    species, rate_consts, inflows, flow_ins, flow_outs = model_export.network_indices(network)

    lines = ["import numpy as np\n"]
    if numba_decoration:
        lines.append("import numba\n")
        lines.append("\n")
        lines.append("@numba.jit(numba.float64[:](numba.float64,numba.float64[:],numba.float64[:]),\n"
                   "\tlocals={'P': numba.float64[:],'F': numba.float64[:,:]},nopython=True)\n")
    lines.append("def model_function(time, S, k):\n")
    lines.append("\n")
    lines.append("\tP = np.zeros(len(S))\n")
    lines.append("\n")

    lines.append("\t")
    lines.append(flow_profile_text)
    lines.append("\n")

    lines.append("\n")
    lines.append("\t")
    lines.append("idx = np.abs(F[0] - time).argmin()")
    lines.append("\n")
    lines.append("\n")
    lines.append("\t")
    lines.append("C = F[:-1,idx]")
    lines.append("\n")
    lines.append("\n")
    lines.append("\t")
    lines.append("k[-1] = F[-1,idx]")
    lines.append("\n")

    lines.append("\n")
    for e in eq_text:
        lines.append("\t{}\n".format(e))

    lines.append("\treturn P\n")
    lines.append("\n")
    lines.append("def wrapper_function(time, S, k):\n")
    lines.append("\treturn model_function(time, S, k)\n")
    lines.append("\n")
    lines.append("species = {")

    for k in species:
        idx = get_index(species[k])
        lines.append("'{}':{},".format(k,idx))
    lines.append("}\n")
    lines.append("\n")

    lines.append("reactions = {")
    for k in rate_consts:
        idx = get_index(rate_consts[k])
        lines.append("'{}':{},".format(k,idx))

    lines.append("}\n")
    lines.append("\n")
    lines.append("inputs = {")
    for k in inflows:
        idx = get_index(inflows[k])
        lines.append("'{}':{},".format(k,idx))

    lines.append("}\n")
    lines.append("\n")
    lines.append("k = np.zeros(max(reactions.values())+1) # rate constants\n")
    lines.append("\n")
    lines.append("S = np.zeros(len(species)) # initial concentrations\n")
    lines.append("\n")

    lines.append("C = np.zeros(len(inputs)) # input concentrations\n")
    lines.append("\n")

    text = ''
    for l in lines:
        text += l

    return text


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
