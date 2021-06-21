import numpy as np
from NorthNet.data_processing import signal_processing
from NorthNet.data_processing import interpolations
from NorthNet import Classes

def combine_datasets(datasets, time_dependent = False):
    '''
    Parameters
    ----------
    datasets: list (len = 2) of NorthNet DataReport objects
        Two datasets with identical conditions to be merged.
    Returns
    -------
    merged_dataset: NorthNet DataReport object.
        Merged datasets with linearly interpolated values
        in new time axis.
    '''
    from NorthNet import Classes

    # check if one data set is void of data,
    # if so, return the other.
    del_idx = False
    for c,d in enumerate(datasets):
        if len(d.series_values) == 0:
            del_idx = c

    if del_idx:
        del datasets[del_idx]
        return datasets[0]

    # create an empty data report
    merged_dataset = Classes.DataReport()
    # Add in information from the first data report in the list
    merged_dataset.experiment_code = datasets[0].experiment_code
    merged_dataset.conditions = datasets[0].conditions
    merged_dataset.analysis_details = datasets[0].analysis_details
    merged_dataset.series_unit = datasets[0].series_unit

    if time_dependent:
        # If the data are time independent, interpolations are
        # performed so that the two datasets share identical time axes.
        # The interpolation is based on the union of the time axes of the
        # two data sets, excluding the values outside of each others'
        # range.
        from scipy import interpolate

        t_a = datasets[0].series_values
        t_b = datasets[1].series_values

        # merge timepoints for the two datasets
        combined_times = np.hstack((t_a, t_b))
        time_axis = np.sort(combined_times)

        # Get the interpolation range
        max_time = min((np.amax(t_a), np.amax(t_b)))
        min_time = max((np.amin(t_a), np.amin(t_b)))
        slice_inds = np.where((time_axis > min_time)&
                              (time_axis < max_time))[0]

        merged_dataset.series_values = time_axis[slice_inds]

        new_data = {}
        for d in datasets:
            for comp in d.data:
                f = interpolate.interp1d(d.series_values,
                                         d.data[comp], kind = "linear")
                new_data[comp] = f(merged_dataset.series_values)

    elif len(datasets[0].series_values) == len(datasets[1].series_values):
        # If the data are considered time-independent and of the same length,
        # measurements between data sets can be considered to be at similar
        # sampling points.
        merged_dataset.series_values = datasets[0].series_values
        new_data = {}

        for d in datasets:
            for comp in d.data:
                new_data[comp] = d.data[comp]
    else:
        # measurements between the two data sets will be combined
        # with the longer one trimmed to the length of the shorter one
        # find the shortest data set
        shortest = min(datasets, key = lambda x:len(x.series_values))
        # Create new data container with the longest dataset trimmed to be
        # the same length as the shortest
        new_data_length = len(shortest.series_values)
        new_data = {}
        for d in datasets:
            for comp in d.data:
                new_data[comp] = d.data[comp][:new_data_length]

        min_val = min(datasets[0].series_values[0],datasets[1].series_values[0])
        merged_dataset.series_values = np.arange(min_val,min_val+new_data_length,1)

    # Put the merged data into the new DataReport object.
    merged_dataset.data = new_data

    return merged_dataset

def get_mass_in_reactor(dataset, flow_profiles, inputs,
                        time_key = 'flow_time/ s'):

    '''
    Parameters
    ----------
    dataset: NorthNet DataReport object
        Data set for mass balance calculation
    flow_profiles: dict
        Dictionary of flow inputs into the experiment.
    flow_keys: list of tuples
        (c,d,e) c = concentration field, d = flow key, e = key to molecular_masses

    Returns
    -------
    mass_in_reactor: numpy 1D array
        Total mass of carbon in the reactor over time
        (same time axis as the input DataReport)
    '''
    import numpy as np
    from NorthNet import info_params

    # get the total flow rate over time from the flow profiles.
    flow_time = flow_profiles[time_key]
    # container for net flow to be added to.
    net_flow = np.zeros(len(flow_profiles[time_key]))
    # iterate through the flow profiles to create an array of the total flow
    # rate for the experiment over time.
    # Further use of this array assumes its values are in seconds.
    for f in flow_profiles:
        if 'flow' in f and not 'time' in f:
            # Conversion to L/s from uL/ h
            net_flow += flow_profiles[f][:len(flow_time)]/(60*60*10e6)

    # get the reactor volume
    # (assumes volume is given in micro litres and is converted to L here)
    reactor_vol = dataset.conditions["reactor_volume/ uL"][0]/10e6

    # Get delta t for mass integration over time.
    deltat = flow_time[1] - flow_time[0]

    # Create empty array for total carbon inputs at each time point.
    mass_input = np.zeros(len(flow_time))
    for i in inputs:
        conc_key = i[0]
        flow_key = i[1]
        # dictionary of molecular masses keyed by SMILES strings
        Mr = info_params.molecular_masses[i[2]]
        # (Mr*flow rate)/(total flow rate)
        mass_input += Mr*flow_profiles[conc_key]*(
                        flow_profiles[flow_key][:len(flow_time)]/
                        (60*60*10e6))*deltat

    # Container for calculated mass in the reactor.
    mass_in_reactor = np.zeros(len(flow_time))
    # Initialise mass in the reactor
    mass_in_reactor[0] = np.nan_to_num(reactor_vol*mass_input[0]/(deltat*net_flow[0]))

    # iterate over inputs and outputs over time to get mass in reactor over time
    for x in range(0,len(mass_in_reactor)-1):
        mass_out = (mass_in_reactor[x]/reactor_vol)*net_flow[x]*deltat
        delta_mass = mass_input[x] - mass_out
        mass_in_reactor[x+1] = mass_in_reactor[x] + delta_mass

    # find where calculated time points are the same as the
    inds = np.where((flow_time > dataset.series_values[0])&(flow_time < dataset.series_values[-1]))[0]
    return mass_in_reactor[inds]

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
