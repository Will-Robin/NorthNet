def interpolate_traces(x,y, length = 1000):
    '''
    For interpolating data
    Parameters
    ----------
    x: numpy array
        X values for the data series to be interpolated.
    y: numpy array
        Y values for the data series to be interpolated.
    length: int
        Number of points to be interpolated

    Returns
    -------

    xnew: numpy array
        Array of the interpolated series' x values.
    ynew: numpy array
        Array of the interpolated series' y values.
    '''

    f = interpolate.interp1d(x,y, kind = "linear")

    xnew = np.linspace(x[0],x[-1], num = length)
    ynew = f(xnew)

    return xnew, ynew
