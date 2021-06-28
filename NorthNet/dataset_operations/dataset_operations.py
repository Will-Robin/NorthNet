import numpy as np

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
