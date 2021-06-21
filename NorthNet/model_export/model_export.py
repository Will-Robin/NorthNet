from NorthNet.model_export import model_export
from NorthNet import info_params
from NorthNet import Classes
import numpy as np


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
