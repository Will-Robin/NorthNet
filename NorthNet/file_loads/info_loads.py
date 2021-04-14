import numpy as np
from NorthNet import Classes

def get_Series_Information(info_file):
    info_dict = {}
    with open(info_file, "r") as f:
        for line in f:
            out = line.strip("\n").split(",")
            info_dict[out[0]] = [x for x in out[1:] if x != ""]

    return info_dict

def import_Experiment_information(fname):
    '''
    Parameters
    ----------
    Get experiment parameters and paths to files
    fname: str
        path to experiment information file
    Returns
    -------
    exp_inf: dict
        Dictionary of Experiment_Information objects. Key: experiment name
    '''

    exp_info = {}
    process_line = lambda x:x.strip("\n").split(",")
    with open(fname) as f:
        for line in f:
            if "Experiment_code" in line:
                header = process_line(line)
            else:
                x = process_line(line)
                exp_info[x[0]] = Classes.Experiment_Information(x[0], x[1], {k:float(v) for k,v in zip(header[3:], x[3:])}, x[2])

    return exp_info

def get_reaction_template_colours(fname, delimiter  = '\t'):

    r_colours = {}
    with open(fname, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = [x for x in line.strip("\n").split(delimiter) if x != '']
                r_colours[ins[0]] = ins[-1]

    return r_colours

def get_generation_protocol(fname):
    '''
    Needs updating
    '''
    with open(fname, 'r') as f:
        for line in f:
            if 'generation_protocol' in line:
                line  = next(f)
                output = line.strip("\n").split(",")[1:]
    output = [x for x in output if x != '']
    return output

def get_environment_matrix(exp_info, input_header,labels, flow_profs):
    from NetFit import dataset_operations as d_ops
    param_keys = [*exp_info[[*exp_info][0]].parameters]
    environments = np.zeros((len(labels),len(input_header)+len(param_keys)))
    for l in range(0,len(labels)):
        ex_k = labels[l].split(";")[0]
        time = float(labels[l].split(";")[1])
        fl_pr = flow_profs[ex_k]

        ins = np.array([exp_info[ex_k].parameters[k] for k in exp_info[ex_k].parameters])
        environments[l] = np.hstack((ins,d_ops.instantaneous_input(time,fl_pr,input_header)))
    return environments
