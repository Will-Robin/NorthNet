import numpy as np

def sum_squares_error(data,x):
    '''
    Calculated sum squared difference between data and x
    data, x: numpy arrays
    '''
    return np.sum((data-x)*(data-x))

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

def polar_to_xy(magnitude, angle):
    '''
    Convert polar coordinates to cartesian
    magnitude: float
    angle: float (in degrees)

    returns: 2-tuple of floats
        x, y
    '''


    x = magnitude*np.cos(np.deg2rad(angle))
    y = magnitude*np.sin(np.deg2rad(angle))

    return x,y

def circle(r, origin):
    '''
    Calculate a circle of radius r and centre origin

    r: float
    origin: 2-tuple of floats

    returns numpy arrays for x and y coordinates
    '''


    theta = np.linspace(0, 2*np.pi, 100)

    x1 = r*np.cos(theta)
    x2 = r*np.sin(theta)

    return x1+origin[0],x2+origin[1]

def rect(x1, y1, x2, y2):
    '''
    get rectangle bounds
    '''
    a = (y1 - y2) / (x1 - x2)
    b = y1 - a * x1
    return (a, b)

def curve_between_points(source,target,centre):
    '''
    create a curve between source and target, moving towards and then away
    from centre

    source,target,centre: 2-tuples of floats
    returns: x_p, y_p list of x and y coordinates for the line
    '''
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
