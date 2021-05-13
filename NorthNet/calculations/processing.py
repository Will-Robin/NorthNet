import numpy as np

def pre_process_cluster_data(data, clusters):

    mod_clusters = {}

    # pre-processing to remove empty entries
    data_mat = np.zeros(len(data))
    for c in clusters:
        data_mat += clusters[c]
    data_mat += data

    idx = np.where(data_mat == 0)[0]

    for c in clusters:
        mod_clusters[c] = np.delete(clusters[c], idx, axis=0)

    mod_data = np.delete(data, idx, axis=0)

    return mod_data, mod_clusters
