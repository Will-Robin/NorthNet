def Euclidean_distances(t_lag_corr_mat):
    '''
    Conversion of time-lagged correlation matrix to a Euclidean distance matrix.
    (see https://en.wikipedia.org/wiki/Euclidean_distance_matrix).
    After Arkin and Ross, J. Phys. Chem. 1995, 99, 970-979.

    Parameters
    ----------
    t_lag_corr_mat: 3D numpy array
        see time_lag_corr_mat() function.

    Returns
    -------
    Euc_dist_mat: 2D numpy array
        A Euclidean distance matrix.

    '''

    absmax_tlcm = np.zeros((len(t_lag_corr_mat),len(t_lag_corr_mat)))
    for d in range(0,len(t_lag_corr_mat)):
        for e in range(0,len(t_lag_corr_mat[d])):
            absmax_tlcm[d,e] = np.amax(np.abs(t_lag_corr_mat[d,e]))

    Euc_dist_mat = np.nan_to_num(np.sqrt(2.0)*np.sqrt(1.0-absmax_tlcm))+1e-10

    return Euc_dist_mat

def multidimensional_scaling(Euc_dist_mat):
    '''
    Multidimensional scaling of a Euclidean distance matrix. Optimises Euclidean
    distances. After Arkin and Ross, J. Phys. Chem. 1995, 99, 970-979.

    Parameters:
    -----------
    Euc_dist_mat: 2D numpy array
        See Euclidean_distances() function.
    '''

    H = np.eye(len(Euc_dist_mat)) - np.ones((len(Euc_dist_mat), len(Euc_dist_mat)))/len(Euc_dist_mat)

    B = -H.dot(Euc_dist_mat**2).dot(H)/2
    eigvals, eigvecs = np.linalg.eigh(B)

    inds   = np.argsort(eigvals)[::-1]
    eigvals = eigvals[inds]
    eigvecs = eigvecs[:,inds]

    # Compute the coordinates using positive-eigenvalued components only
    w, = np.where(eigvals > 0)
    L  = np.diag(np.sqrt(eigvals[w]))
    V  = eigvecs[:,w]
    Y  = V.dot(L)
    k = 1
    a1k = np.sum(eigvals[:k])/np.sum(np.abs(eigvals[:k]))
    a2k = np.sum(eigvals[:k]**2)/np.sum(eigvals[:k]**2)

    return np.nan_to_num(Y)
