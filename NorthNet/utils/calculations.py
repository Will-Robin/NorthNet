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
