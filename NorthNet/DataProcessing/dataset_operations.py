import numpy as np
from NorthNet.DataProcessing import signal_processing
from NorthNet.DataProcessing import interpolations
from NorthNet import Classes

def get_mass_in_reactor(dataset, flow_profiles, inputs):

    '''
    flow_keys: list of tuples
        (c,d,e) c = concentration field, d = flow key, e = key to molecular_masses
    '''
    time = flow_profiles['flow_time/ s']
    net_flow = np.zeros(len(flow_profiles['flow_time/ s']))
    for f in flow_profiles:
        if 'flow' in f and not 'time' in f:
            net_flow += flow_profiles[f][:len(time)]/(60*60*10e6)

    reactor_vol = dataset.conditions["reactor_volume/ uL"][0]/10e6

    deltat = time[1] - time[0]

    mass_input = np.zeros(len(time))
    for i in inputs:
        for f in flow_profiles:
            if i in f and 'flow' in f:
                flow_k = f
            elif i in f and not 'flow' in f:
                inp_k = f

        mass_input += info_params.molecular_masses[i]*flow_profiles[inp_k]*(flow_profiles[flow_k][:len(time)]/(60*60*10e6))*deltat

    mass_in_reactor = np.zeros(len(time))

    mass_in_reactor[0] = np.nan_to_num(reactor_vol*mass_input[0]/(deltat*net_flow[0]))

    for x in range(0,len(mass_in_reactor)-1):
        mass_out = (mass_in_reactor[x]/reactor_vol)*net_flow[x]*deltat
        delta_mass = mass_input[x] - mass_out
        mass_in_reactor[x+1] = mass_in_reactor[x] + delta_mass

    inds = np.where((time > dataset.time[0])&(time < dataset.time[-1]))[0]
    return time[inds], mass_in_reactor[inds]

def get_amplitudes(dataset, write_f = False):
    '''
    Get amplitudes of the variables in a Dataset objects via Fourier transform.

    Parameters
    ----------
    dataset: Dataset object
        Object containig data.

    Returns
    -------
    amps_dict: dict
        Dictionary of amplitudes of species indexed by their names from Dataset.
    '''

    deps = [*dataset.dependents]
    time = dataset.time
    amps_dict = {}

    for ds in deps:
        amps = signal_processing.fourier_amplitude(time,dataset.dependents[ds])
        amps_dict[ds] = amps

    if write_f:
        with open("{}_amplitudes.csv".format(dataset.name), "w") as f:
            f.write("{},{}\n".format("species","amplitude/ M"))
            for a in amps_dict:
                f.write("{},{}\n".format(a,amps_dict[a]))

    return amps_dict

def get_amplitudes_as_percent(dataset, write_f = False):
    '''
    Get amplitudes of the variables in a Dataset objects via Fourier transform.

    Parameters
    ----------
    dataset: Dataset object
        Object containig data.

    Returns
    -------
    amps_dict: dict
        Dictionary of amplitudes of species indexed by their names from Dataset.
    '''

    deps = [*dataset.dependents]
    time = dataset.time
    amps_dict = {}

    for ds in deps:
        amps = signal_processing.fourier_amplitude(time,dataset.dependents[ds])
        amps_dict[ds] = 100*amps/np.average(dataset.dependents[ds])

    if write_f:
        with open("{}_amplitudes.csv".format(dataset.name), "w") as f:
            f.write("{},{}\n".format("species","amplitude/ M"))
            for a in amps_dict:
                f.write("{},{}\n".format(a,amps_dict[a]))

    return amps_dict

def get_wave_parameter_by_fitting(dataset, key = "amplitude", write_f = False):
    '''
    Get amplitudes of the variables in a Dataset objects via Fourier transform.

    Parameters
    ----------
    dataset: Dataset object
        Object containig data.

    key: str
        Choose from "amplitude","phase","offset","period".

    write_f: bool
        Whether to write an output file or not.

    Returns
    -------
    amps_dict: dict
        Dictionary of amplitudes of species indexed by their names from Dataset.
    '''

    deps = [*dataset.dependents]
    time = dataset.time
    output_dict = {}

    for ds in deps:
        params = fitting_functions.fit_sine_wave(dataset.time, dataset.dependents[ds])
        output_dict[ds] = params[key]

    if write_f:
        with open("{}_{}.csv".format(dataset.name, key), "w") as f:
            f.write("{},{}\n".format("species","{}/ SI_unit".format(key)))
            for a in output_dict:
                f.write("{},{}\n".format(a,output_dict[a]))

    return output_dict

def data_to_amplitude_correlations(multiple_data_sets,info_dict):
    '''
    Get amplitudes of dependents from a collection of datasets and combine the
    information.

    Parameters
    ----------
    multiple_data_sets: dict
        dictionary containing Dataset objects
    info_dict: dict
        Dictionary containing information (drive frequency and residence times)
        for datasets.

    Returns
    -------
    amps_dict: dict
        Dictionary of dictionaries amplitudes indexed by dataset.

    '''
    amps_dict = {}
    dr_fr = []
    res_tim = []
    for c,d in enumerate(multiple_data_sets):

        deps = [*multiple_data_sets[d].dependents]

        for ds in deps:
            amps = signal_processing.fourier_amplitude(multiple_data_sets[d].time,multiple_data_sets[d].dependents[ds])
            nw_dr_fr = info_dict[d][1]
            nw_res_time  = info_dict[d][0]
            if ds in amps_dict:
                amps_dict[ds]["amps"].append(amps)
                amps_dict[ds]["drive_freq"].append(nw_dr_fr)
                amps_dict[ds]["residence_time"].append(nw_res_time)
            else:
                amps_dict[ds] = {}
                amps_dict[ds]["amps"] = []
                amps_dict[ds]["drive_freq"] = []
                amps_dict[ds]["residence_time"] = []

                amps_dict[ds]["amps"].append(amps)
                amps_dict[ds]["drive_freq"].append(nw_dr_fr)
                amps_dict[ds]["residence_time"].append(nw_res_time)

    return amps_dict

def get_input_profile(data_set):
    return data_set.get_input_profile()

def combine_datasets(data_sets, conditions, interpolation_length = 100, time_independent = False, x_axis_key = 'time/ s'):

    '''
    Combines two separate datasets into one dataset object

    Parameters
    ---------
    data_sets: list
        List of data sets to be combined. All datasets must have similarly
        spaced x-axis values.
    conditions: dict
        Dictionary of conditions to be input into the new dataset.

    Returns
    -------
        Dataset object combinging the two datasets named after the last dataset
        in the list with specified conditions.
    '''
    temp_interp = {}

    if not time_independent:
        t_max = 1e100
        t_min = 0

        for x in data_sets:

            if x.time[-1] < t_max:
                t_max = x.time[-1]
            if x.time[0] > t_min:
                t_min = x.time[0]

        for x in data_sets:

            deps = sorted([*x.dependents])

            for d in deps:
                inds = np.where((x.time >= t_min)&(x.time <= t_max))[0]
                if len(inds) == 0:
                    pass
                else:
                    new_x, new_y = interpolations.interpolate_traces(x.time[inds],x.dependents[d][inds], length = interpolation_length)
                    temp_interp[d] = new_y
        else:
            for x in data_sets:
                deps = sorted([*x.dependents])
                for d in deps:
                    new_x, new_y = interpolations.interpolate_traces(x.time,x.dependents[d], length = interpolation_length)
                    temp_interp[d] = new_y

        temp_interp[x_axis_key] = np.linspace(t_min,t_max, num = interpolation_length)

        return Classes.Dataset("{}_interpolated".format(x.name), conditions, temp_interp)

    else:
        deps = {}
        for x in data_sets:
            deps[x_axis_key] = x.time
            for d in x.dependents:
                deps[d] = x.dependents[d]
        print('combining', x_axis_key)
        return Classes.Dataset("{}".format(x.name), conditions, deps, time_key = x_axis_key, get_flow_profile = False)

def add_flow_segment(profile_time, profile, dataset):
    profile_time = np.array(profile_time)
    profile = np.array(profile)

    clip = np.where((profile_time>dataset.time[0]) & (profile_time < dataset.time[-1]))

    new_x, new_y = interpolations.interpolate_traces(profile_time[clip],profile[clip], length = 1000)

    setattr(dataset, "drive", new_y)

def split_data_into_sections(dataset, time_split = 1000):
    '''
    A function to split a dataset up into sections based on sections of the
    data being separated by the time_split value.

    Parameters
    ----------
    dataset: NetFit Dataset object
        Dataset to be split into sections.
    time_split: float
        Time separating sections of dataset.

    Returns
    -------
    store: list of dataset objects
        The dataset split up into components as a list.
    '''

    store = []

    dif  = np.diff(dataset.time)
    inds = np.where(dif > time_split)[0]
    time_brackets = dataset.time[inds]

    time_brackets = np.hstack((0,time_brackets))
    time_brackets = np.hstack((time_brackets, dataset.time[-1]))

    for i in range(0,len(time_brackets)-1):

        data_insert = {}
        inds = np.where((dataset.time > time_brackets[i])&(dataset.time <= time_brackets[i+1]))[0]
        data_insert["time/ s"] = dataset.time[inds]

        for d in dataset.dependents:
            data_insert[d]  = dataset.dependents[d][inds]

        new_dset = Classes.Dataset(dataset.name + "_{}".format(i), dataset.conditions, data_insert)

        store.append(new_dset)

    return store

def instantaneous_input(time,flow_profile,input_header):

    t_pos =  time - time%2 # time at datapoint was taken, rounded to lowest 2 because flow profiles are in 2 s time intervals

    inp_time = np.where(flow_profile["flow_time/ s"] == int(t_pos))[0] # index of flows at t_pos

    inp_arr = np.zeros(len(input_header))
    tot_flow = np.zeros(len(flow_profile["flow_time/ s"]))
    inst_cond = np.zeros(len(input_header))

    for f in flow_profile:
        if f == "flow_time/ s":
            pass
        elif "/ M" in f and not "water" in f:
            pos = input_header.index(f.split("/")[0])
            inp_arr[pos] = flow_profile[f][0]
        elif "/h" in f:
            tot_flow += flow_profile[f][inp_time]

    for f in flow_profile:
        if "/h" in f and f.split("_")[0] in input_header:
            pos = input_header.index(f.split("_")[0])
            inst_cond[pos] = (inp_arr[pos]* flow_profile[f][inp_time])/ tot_flow[inp_time]

    return inst_cond
