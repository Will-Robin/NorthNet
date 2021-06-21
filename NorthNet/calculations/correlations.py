import numpy as np

def autocorrelation(time, trace):
    '''
    Autocorrelation for a time trace.
    '''

    roll_vals = np.arange(-(len(time))+1,len(time))
    t_step = time[1]-time[0]

    t_lags = np.zeros(len(roll_vals))

    for r in range(0,len(roll_vals)):
        t_lags[r] = roll_vals[r]*t_step


    x_i = trace - trace.mean()
    con = np.correlate(x_i,x_i, mode = "full")
    ccor = con/(len(time)*x_i.std()*x_i.std())

    return t_lags, ccor

def correlation_matrix(arr):
    '''
    axis 0 is dataset

    '''

    corr_mat = np.zeros((len(arr),len(arr)))

    cent = np.zeros(arr.shape)
    for x in range(0,len(arr)):
        cent[x] = arr[x] - np.mean(arr[x])

    for x in range(0,len(cent)):
        x_i = cent[x]
        for y in range(0,len(cent)):
            x_j = cent[y]
            con = np.correlate(x_i,x_j, mode = "valid")
            corr_mat[x,y] = np.nan_to_num(con/(len(x_i)*np.std(x_i, ddof = 1)*np.std(x_j, ddof = 1)))

    return corr_mat

def time_lag_corr_mat(dataset):
    '''
    Create a time-lagged correlation matrix from a dataset object.

    Parameters:
    -----------
    dataset: Dataset object
        Object containing the data. Note that the data mist be equally spaced in
        the time dimension.

    Returns
    -------
    t_lag_corr_mat: 3D numpy array
        3D numpy matrix containing cross correlations of species.

        axis 0: length of dependents in Dataset object, indexed by
        [*dataset.depedents].
        axis 1: length of dependents in Dataset object, indexed by
        [*dataset.depedents].
        axis 2: length of time values. Indexes by time-lag values
        (see t_lags variable)

    '''

    time = dataset.time
    signals = np.array([])

    deps = [*dataset.dependents]

    roll_vals = np.arange(-(len(time))+1,len(time))
    t_step = time[1]-time[0]

    t_lags = np.zeros(len(roll_vals))

    for r in range(0,len(roll_vals)):
        t_lags[r] = roll_vals[r]*t_step

    t_lag_corr_mat = np.zeros((len(deps),len(deps),len(t_lags)))

    for d in range(0,len(deps)):

        x_i = dataset.dependents[deps[d]] - dataset.dependents[deps[d]].mean()

        for e in range(0,len(deps)):
            x_j = dataset.dependents[deps[e]] - dataset.dependents[deps[e]].mean()

            con = np.correlate(x_i,x_j, mode = "full")
            ccor = con/(len(time)*x_i.std()*x_j.std())

            t_lag_corr_mat[d,e] = np.nan_to_num(ccor)

    return t_lag_corr_mat

def t_lag_corr_to_drive(dataset, modulated_input_key = "dihydroxyacetone_concentration/ M",
                        flow_profile_key = "dihydroxyacetone_flow_profile/ uL/h"):

    if modulated_input_key not in dataset.conditions:
        print([*dataset.conditions])
    elif flow_profile_key not in dataset.input_flows:
        print([*dataset.input_flows])

    time = dataset.time
    deps = [*dataset.dependents]
    modulated_input = dataset.conditions[modulated_input_key]*dataset.input_flows[flow_profile_key]/dataset.net_flow

    roll_vals = np.arange(-(len(time))+1,len(time))
    t_step = time[1]-time[0]

    t_lag_corr_mat = np.zeros((len(deps),len(roll_vals)))
    x_i = modulated_input - modulated_input.mean()

    for e in range(0,len(deps)):
        x_j = dataset.dependents[deps[e]] - dataset.dependents[deps[e]].mean()

        con = np.correlate(x_i,x_j, mode = "full")
        ccor = con/(len(time)*x_i.std()*x_j.std())

        t_lag_corr_mat[e] = np.nan_to_num(ccor)

    return t_lag_corr_mat


def time_lags_to_drive(t_lag_corr_mat,time):

    roll_vals = np.arange(-(len(time))+1,len(time))
    t_step = time[1]-time[0]

    t_lags = np.zeros(len(roll_vals))

    for r in range(0,len(roll_vals)):
        t_lags[r] = roll_vals[r]*t_step

    lag_mat = np.zeros((len(t_lag_corr_mat)))

    for d in range(0,len(t_lag_corr_mat)):
        test_line = np.round(t_lag_corr_mat[d],3)
        abs_line = np.abs(test_line)
        ind = np.where((abs_line == np.amax(abs_line)))[0]
        lag_mat[d]  = t_lags[ind[0]]

    return lag_mat

def time_lags(t_lag_corr_mat,time):
    '''
    Get time lags from time-lagged correlation matrix.

    Parameters:
    -----------
        t_lag_corr_mat: 3D numpy array
            see time_lag_corr_mat() function.
        time: numpy array
            array of time values used to create t_lag_corr_mat.

    Returns
    -------
    lag_mat: 2D numpy array
        Array of time lag values indexed to original data in the same way as
        t_lag_corr_mat
    '''

    roll_vals = np.arange(-(len(time))+1,len(time))
    t_step = time[1]-time[0]

    t_lags = np.zeros(len(roll_vals))

    for r in range(0,len(roll_vals)):
        t_lags[r] = roll_vals[r]*t_step

    lag_mat = np.zeros((len(t_lag_corr_mat),len(t_lag_corr_mat)))
    to_pass = []
    for d in range(0,len(t_lag_corr_mat)):
        to_pass.append(d)
        for e in range(0,len(t_lag_corr_mat[d])):
            if e in to_pass:
                pass
            else:
                test_line = np.round(t_lag_corr_mat[d,e],3)
                abs_line = np.abs(test_line)
                ind = np.where((abs_line == np.amax(abs_line)))[0]

                lag_mat[d,e]  = t_lags[ind[0]]

    return lag_mat

def t_lag_corr_mat_to_correlation_matrix(t_lag_corr_mat):
    '''
    Create correlation matrix from a time-lagged correlation matrix.
    Parameters
    ----------
    t_lag_corr_mat: 3D numpy array
        see time_lag_corr_mat() function.

    Returns
    -------
    corr_mat: 2D numpy array
        Array of correlation coefficient values indexed to original data in the
        same way as t_lag_corr_mat
    '''

    corr_mat = np.zeros((len(t_lag_corr_mat),len(t_lag_corr_mat[0])))
    to_pass = []
    for d in range(0,len(t_lag_corr_mat)):
        to_pass.append(d)
        for e in range(0,len(t_lag_corr_mat[d])):
            if e in to_pass:
                pass
            else:
                test_line = np.round(t_lag_corr_mat[d,e],3)
                abs_line = np.abs(test_line)
                ind = np.where((abs_line == np.amax(abs_line)))[0]
                corr_mat[d,e]  = t_lag_corr_mat[d,e,ind[0]]

    return np.nan_to_num(corr_mat)

def connection_table(t_lag_corr_mat, deps):
    '''Create connection table'''
    groups = {}

    cor_map = np.zeros((len(t_lag_corr_mat),len(t_lag_corr_mat)))
    for d in range(0,len(t_lag_corr_mat)):
        for e in range(0,len(t_lag_corr_mat[d])):

            max_corr = np.amax(t_lag_corr_mat[d,e])
            min_corr = np.amin(t_lag_corr_mat[d,e])

            if np.abs(min_corr) > max_corr:
                strong_corr = min_corr
            else:
                strong_corr = max_corr

            if abs(strong_corr) > 0.1:

                gr = sorted([deps[d],deps[e]])
                if "{},{}".format(gr[0],gr[1]) in groups:
                    if strong_corr > groups["{},{}".format(gr[0],gr[1])]:
                        groups["{},{}".format(gr[0],gr[1])] = strong_corr
                    else:
                        pass
                else:
                    groups["{},{}".format(gr[0],gr[1])] = strong_corr
            else:
                pass

    return groups
