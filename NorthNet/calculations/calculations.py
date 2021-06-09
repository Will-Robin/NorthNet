def calculate_weighted_vector(clusters, weights):
    import numpy as np
    sum = np.zeros(len(clusters[[*clusters][0]]))
    for c,clust in enumerate(clusters):
        sum += weights[c]*clusters[clust]
    return sum

def _error(data,x):
    import numpy as np
    return np.sum((data-x)*(data-x))

def compare(weights, data, clusters):
    from NorthNet.calculations import calculations
    x = calculations.calculate_weighted_vector(clusters, weights)
    return _error(data,x)

def test_vector(newvec, data, reconst):
    x = newvec + reconst
    return _error(data,x)

def sine_wave(time, period, amplitude, phase, offset):
    '''
    For generating a sine wave flow profile.

    Parameters
    ----------
    period: float
        The period of the sine wave in seconds.

    amplitude: float
        The amplitude of the sine wave (any unit: think of output).

    phase:
        Phase shift

    offset: float
        The centre of the sine wave.

    Returns
    -------
    wave: 1Darray
        An array of flow rate values.
    '''

    wave = amplitude*np.sin(2*np.pi*time/period + phase) + offset

    return wave

def compare_model_data(k, data, S, model):
    '''
    For determining the error between a calculated model and data. Kind of old
    and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    time,calc = integrate_model(data.time, data.initial, k, model)

    spec_inds = {v:k for k,v in zip([*model.species],model.species.values())}

    error = np.zeros([len(calc[:,0])])

    for c in range(0,len(calc[0])):
        if spec_inds[c]+" M" in [*data.dependents]:
            error += error + (calc[:,c] - data.dependents[spec_inds[c]+" M"])*(calc[:,c] - data.dependents[spec_inds[c]+" M"])/np.amax(data.dependents[spec_inds[c]+" M"])

    return np.sum(error)

def error_calculation(data, calc):
    ''' Method for calculating the error between the data and calculated fits.
    Kind of old and not used in most things this module is used for.
    Should be modified and documented later.
    '''

    err = ((data - calc)/data)*((data - calc)/data)

    error = np.sum(err)/len(data)

    return error

def polar_to_xy(magnitude, angle):
    import numpy as np

    x = magnitude*np.cos(np.deg2rad(angle))
    y = magnitude*np.sin(np.deg2rad(angle))

    return x,y

def circle(r, origin):

    theta = np.linspace(0, 2*np.pi, 100)

    x1 = r*np.cos(theta)
    x2 = r*np.sin(theta)

    return x1+origin[0],x2+origin[1]

def rect(x1, y1, x2, y2):
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a * x1
    return (a, b)

def curve_between_points(source,target,centre):
    from NorthNet.calculations import calculations

    (x1, y1, x2, y2) = (source[0],source[1],centre[0],centre[1])
    (a1, b1) = calculations.rect(source[0],source[1], centre[0],centre[1])
    (a2, b2) = calculations.rect(target[0],target[1], centre[0],centre[1])
    x_p = []
    y_p = []
    for i in range(0, 1000):
        if x1 == x2:
            continue
        else:
            (a, b) = calculations.rect(x1, y1, x2, y2)
        x = i*(x2 - x1)/1000 + x1
        y = a*x + b
        x_p.append(x)
        y_p.append(y)
        x1 += (centre[0] - source[0])/1000
        y1 = a1*x1 + b1
        x2 += (target[0] - centre[0])/1000
        y2 = a2*x2 + b2

    return x_p, y_p
