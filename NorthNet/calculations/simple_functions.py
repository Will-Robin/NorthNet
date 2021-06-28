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

def get_rt_from_header(element):
    from NorthNet.misc import simple_functions as s_f

    if s_f.isfloat(element):
        position = float(element)
    else:
        spl = element.split('(')[-1].strip(')')
        position = float(spl[0:-1])
    return position
