from NorthNet.calculations import data_fitting
from NorthNet.calculations import processing
from NorthNet.calculations import calculations
import numpy as np

def fit_clusters(data, clusters):
    from scipy.optimize import minimize
    ans = minimize(calculations.compare, np.random.rand(len(clusters)), args = (data,clusters), bounds = [(0,10) for x in range(0,len(clusters))], tol = 1e-9)
    return ans.x

def fit_vectors_to_data(data, clusters, header):

    mod_data, mod_clusters = processing.pre_process_cluster_data(data, clusters)

    weights = data_fitting.fit_clusters(mod_data, mod_clusters)

    x = calculations.calculate_weighted_vector(clusters, weights)
    resid = data-x
    pc_resid = np.abs(np.sum(resid))/np.sum(data)

    return weights, pc_resid

def fit_discrepancy_vector(data, clusters, weights):

    from scipy.optimize import minimize
    reconst = calculations.calculate_weighted_vector(clusters, weights)
    ans = minimize(data_fitting.test_vector, np.random.rand(len(clusters[[*clusters][0]])), args = (data,reconst), bounds = [(0,10) for x in range(0,len(clusters[[*clusters][0]]))], tol = 1e-9)

    return ans.x

def fit_clusters_all_data(data_sets, clusters, header, remove_inds):
    wghts = []
    residuals = []
    tuning_vectors = []
    resid_vecs = []
    lb_dat = []
    for d in data_sets:
        vector_dict = data_fitting.convert_dataset_to_vector(data_sets[d], header)
        for a in vector_dict:
            for i in remove_inds:
                vector_dict[a][i] = 0.0
        for c,a in enumerate(vector_dict):
            phenotype,resid = data_fitting.fit_vectors_to_data(vector_dict[a], clusters, header)
            residuals.append(resid)
            wghts.append(phenotype)
            tuning_vectors.append(data_fitting.fit_discrepancy_vector(vector_dict[a], clusters, phenotype))
            resid_vecs.append(vector_dict[a] - calculations.calculate_weighted_vector(clusters, phenotype))
            lb_dat.append("{};{}".format(d,data_sets[d].time[c]))

    return wghts, residuals, tuning_vectors, resid_vecs, lb_dat

def cross_fit_clusters(clusters):
    wghts = []
    residuals = []
    tuning_vectors = []
    resid_vecs = []
    lb_dat = []
    cross_mat = np.zeros((len(clusters),len(clusters)))
    for c1,d in enumerate(clusters):
        loc_clust = {k:clusters[k] for k in clusters}
        loc_clust[d] = np.ones(len(clusters[d]))
        phenotype,resid = data_fitting.fit_vectors_to_data(clusters[d], loc_clust, header)
        cross_mat[c1] = phenotype
    return cross_mat

def fit_sine_wave(time,signal):
    '''
    For fitting a sine wave to a signal trace

    Paramters
    ---------
    time: 1D numpy array
        Time array
    signal: 1D numpy array
        Intensity of signal over time.

    Returns
    -------
    popt: list
        list of fitting parameters: [period, amplitude, phase, offset]

    '''

    if np.average(signal) == 0.0:
        return {"period":1, "amplitude":0, "phase":0, "offset":0}
    else:
        initial_guess = [(time[-1]-time[0])/3, (np.amax(signal)-np.amin(signal))/2, 0, np.average(signal)]
        lowerbounds   = [0, (np.average(signal)- np.amax(signal))*3, -np.pi, 0.0]
        upperbounds   = [(time[-1]-time[0]), (np.amax(signal)-np.average(signal))*3, np.pi, 2*(np.average(signal))]
        boundarr  = [lowerbounds,upperbounds]

        try:
            popt, pcov = curve_fit(fitting_functions.sine_wave, time, signal, p0=initial_guess, bounds = boundarr)

            if np.average(signal) - popt[1] < 0:
                popt[1] = 0

            return {"period":popt[0], "amplitude":popt[1], "phase":popt[2], "offset":popt[3]}

        except:

            return {"period":1, "amplitude":0, "phase":0, "offset":np.average(signal)}
