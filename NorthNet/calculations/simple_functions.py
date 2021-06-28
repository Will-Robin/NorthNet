def cluster(values, bound = 0.1):

    import numpy as np

    values = np.sort(values)

    cluster = []
    for m in range(0,len(values)):
        if len(cluster) == 0:
            pass
        else:
            clust_av = np.average(cluster)

            if  abs(values[m]-clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(values[m])

    yield cluster

def isfloat(thing):
    '''
    Test if an thing (e.g. str) can be converted to a float.
    '''
    try:
        float(thing)
        return True
    except ValueError:
        return False
