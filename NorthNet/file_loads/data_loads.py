import numpy as np
from NorthNet.information import info_params

def get_data(file, x_axis_key = "time/ s", flow_data = True):
    '''
    From data_fitting_functions.py
    For creating a Dataset object from a report format file.

    Paramters
    ------
    file: str
         Name of the data file in which the data are stored. Must be a .csv file,
         comma delimited and in the correct data repot format.

    Returns
    -------
    Dataset object:
        A dataset object created from the file.
    '''

    condset = []
    readstate = False
    with open(file, "r", encoding = 'latin-1') as f:
        for c,line in enumerate(f):
            if 'Dataset' in line:
                dset_name = line.strip("\n").split(",")[1]
            if "start_conditions" in line:
                readstate = True
                line = next(f)
            if "end_conditions" in line:
                readstate = False
            if readstate:
                newline = line.strip("\n")
                condset.append([x for x in newline.split(",") if x != ""])
    c_out = {}
    for c in condset:
        try:
            c_out[c[0]] = [float(x) for x in c[1:]]
        except:
            c_out[c[0]] = [x for x in c[1:]]

    dataset = []
    with open(file, 'r',encoding='latin1') as f:
        for n,line in enumerate(f,0):
            if "start_data" in line:
                readstate = True
                line = next(f)
            if "end_data" in line:
                readstate = False
            if readstate:
                newline = line.strip("\n")
                dataset.append([x for x in newline.split(",") if x != ""])

    e = [list(i) for i in zip(*dataset)]
    d_out = {}
    for s in e:
        d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

    return Classes.Dataset(dset_name,c_out,d_out, time_key = x_axis_key,  get_flow_profile = flow_data)

def load_data_set(path):
    from NorthNet import dataset_operations as d_ops
    home = os.getcwd()
    t_inep = False
    os.chdir(path)
    lst = os.listdir()
    d_sets = []
    for f in lst:
        if f.endswith(".csv"):
            try:
                dset = f_io.get_data(f)
                d_sets.append(dset)
                print('loaded flow')
            except:
                t_inep = True
                dset = f_io.get_data(f, x_axis_key = "sample number/ n", flow_data = False)
                d_sets.append(dset)
                print('not flow')

    os.chdir(home)
    return d_ops.combine_datasets(d_sets, d_sets[0].conditions, time_independent = t_inep,
                                  interpolation_length = len(d_sets[0].time),
                                  x_axis_key = d_sets[0].independent_name)

def data_from_directory(directory):
    '''
    For importing multiple datasets from a directory into a dictionary indexed
    by file names.

    Parameters
    ----------
    directory: str
        Path to the directory holding the data files.

    Returns
    -------
    data_dict: dict
        A dictionary of dataset objects indexed by the first 8 characters of
        their filenames
    '''

    cwd = os.getcwd()

    os.chdir(directory)

    f_list = [x for x in os.listdir() if x.endswith(".csv")]
    data_dict = {}
    for f in f_list:
        if f[:8] in data_dict:
            data_dict[f[:8]].append(get_data(f))
        else:
            data_dict[f[:8]] = []
            data_dict[f[:8]].append(get_data(f))

    os.chdir(cwd)

    return data_dict

def load_flow_profiles(fname):
    flow_info = {}
    proc_line = lambda x:[z for z in x.strip("\n").split(",") if z != ""]
    with open(fname, "r") as f:
        for line in f:
            ins = proc_line(line)
            if ins[0] in flow_info:
                flow_info[ins[0]][ins[1]] = np.array([float(z) for z in ins[2:]])
            else:
                flow_info[ins[0]] = {}
                flow_info[ins[0]][ins[1]] = np.array([float(z) for z in ins[2:]])
    return flow_info

def load_exp_compound_file(fname, header):
    output = {}
    with open(fname, 'r') as f:
        for c,line in enumerate(f):
            if c == 0:
                ins  = line.strip('/n').split(',')
                f_head = ins[1:]
            else:
                ins = line.strip('/n').split(',')
                data = ins[1:]
                fill_line = [0.0]*len(header)
                for x in range(0,len(f_head)):

                    print(header)
                    print(f_head[x])
                    
                    idx = header.index(f_head[x])
                    fill_line[idx] = data[x]

                output[ins[0]] = fill_line
    return output

def load_cluster_locations(fname):
    output = {}
    proc_line = lambda x:[z for z in x.strip("\n").split(",") if z != ""]
    with open(fname, "r") as f:
        for line in f:
            ins = proc_line(line)
            output[ins[0]] = ins[1:]
    return output

def get_series_vectors(exp_info, sets = [], param = 'offsets',second_folder_path= "Wave_parameters"):
    '''
    Parameters
    ----------
    exp_info: dict
        Dictionary of Experiment_Information objects, key in experiment name
    sets: list
        String keys to exp_info. If empty, all sets loaded
    second_folder_path:
        Folder containing experiment results within main experiment folder.

    Returns
    -------
    output: dict
        Dictionary of compositional vectors: keys: experiment names
        Indices of compounds are the same as in header
    header: list
        List of compounds: indices of vectors in output are same as in header.
    '''

    home = os.getcwd()
    # get composition vectors
    if len(sets) == 0:
        sets = [*exp_info]

    vector_dict = {x+" M":[0.0]*len(sets) for x in info_params.smiles_to_names}

    for count, e in enumerate(sets):
        f_path = Path(exp_info[e].path) / second_folder_path
        os.chdir(f_path)
        f_list = os.listdir()
        for f in f_list:
            if exp_info[e].name[:-1] in f:
                os.chdir(f)
                f2_level = os.listdir()
                for f2 in f2_level:
                    if f2.endswith('.csv'):
                        amps, phases, offsets = f_io.read_wave_parameters_file(f2)
                break

        if param == 'amps':
            for c,a in enumerate(amps):
                vector_dict[a][count] = amps[a]
        if param == 'phases':
            for c,a in enumerate(phases):
                vector_dict[a][count] = phases[a]
        if param == 'offsets':
            for c,a in enumerate(offsets):
                vector_dict[a][count] = offsets[a]

    os.chdir(home)

    stack = np.zeros((len(vector_dict),len(sets)))

    for c,v in enumerate(vector_dict):
        stack[c] = vector_dict[v]

    matrix = stack.T

    output = {}
    for x in range(0,len(matrix)):
        output[sets[x]] = matrix[x]

    header = [*vector_dict]

    return output, header

def get_composition_vectors(exp_info, sets = [], second_folder_path= "Wave_parameters"):
    '''
    Parameters
    ----------
    exp_info: dict
        Dictionary of Experiment_Information objects, key in experiment name
    sets: list
        String keys to exp_info. If empty, all sets loaded
    second_folder_path:
        Folder containing experiment results within main experiment folder.

    Returns
    -------
    output: dict
        Dictionary of compositional vectors: keys: experiment names
        Indices of compounds are the same as in header
    header: list
        List of compounds: indices of vectors in output are same as in header.
    '''

    home = os.getcwd()
    # get composition vectors
    if len(sets) == 0:
        sets = [*exp_info]

    vector_dict = {x+" M":[0.0]*len(sets) for x in info_params.smiles_to_names}

    for count, e in enumerate(sets):
        f_path = Path(exp_info[e].path) / second_folder_path
        os.chdir(f_path)
        f_list = os.listdir()
        for f in f_list:
            if exp_info[e].name[:-1] in f:
                os.chdir(f)
                f2_level = os.listdir()
                for f2 in f2_level:
                    if f2.endswith('.csv'):
                        amps, phases, offsets = f_io.read_wave_parameters_file(f2)
                break

        for c,a in enumerate(offsets):
            vector_dict[a][count] = offsets[a]

    os.chdir(home)

    stack = np.zeros((len(vector_dict),len(sets)))

    for c,v in enumerate(vector_dict):
        stack[c] = vector_dict[v]

    matrix = stack.T

    output = {}
    for x in range(0,len(matrix)):
        output[sets[x]] = matrix[x]

    header = [*vector_dict]

    return output, header

def get_composition_base_vectors(fname):
    clusters = {}
    with open(fname,"r") as f:
        for line in f:
            if "cluster" in line:
                nom = line.strip("\n").split(",")[1]
            if "vector" in line:
                ins = line.strip("\n").split(",")[1:]
                clusters[nom] = np.array([float(x) for x in ins if x != ""])
    return clusters

def read_mass_balance_file(filename,header):
    mass_balances = {}
    labels = {}
    with open(filename, 'r') as f:
        for line in f:
            ins = line.strip('\n').split(',')
            if '_labels' in ins[0]:
                labels[ins[0].split('_')[0]] = [x for x in ins[1:] if x != '']
            if '_percentage_mass' in ins[0]:
                mass_balances[ins[0].split('_')[0]] = [float(x) for x in ins[1:] if x != '']

    output = {m:[] for m in mass_balances}
    for l in labels:
        vec = np.zeros(len(header))
        for c,v in enumerate(labels[l]):
            ind = header.index(v+' M')
            vec[ind] = mass_balances[l][c]

        output[l] = vec

    return output

def get_reaction_expressions(fname):
    rxn_clusters = {}
    process_line = lambda x: x.strip("\n").split(",")
    with open(fname, "r") as f:
        for line in f:
            if 'Chemotype' in line:
                name = line.strip('\n')
            if 'Reaction_class' in line:
                ins = process_line(line)
                header = [x for x in ins[1:] if x != ""]
            if 'Count' in line:
                ins = process_line(line)
                rxn_clusters[name] = np.array([int(x) for x in ins[1:] if x != ""])

    return header,rxn_clusters

def read_wave_parameters_file(fname):
    amps = {}
    phases = {}
    offsets = {}
    with open(fname, "r") as f:
        next(f)
        for line in f:

            x = line.strip("\n").split(",")

            amps[x[0]] = float(x[1])
            phases[x[0]] = float(x[2])
            offsets[x[0]] = float(x[3])

    return amps, phases, offsets

def get_wave_parameters(path):
    cwd = os.getcwd()
    pt = path/"Wave_parameters"
    os.chdir(pt)
    f_list = os.listdir()
    for f in f_list:
        if e.name[:-1] in f:
            os.chdir(f)
            f2_level = os.listdir()
            for f2 in f2_level:
                if f2.endswith('.csv'):
                    amps, phases, offsets = f_io.read_wave_parameters_file(f2)
            break

    os.chdir(cwd)

    return amps, phases, offsets
