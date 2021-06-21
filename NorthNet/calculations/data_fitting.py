import numpy as np

def fit_sine_wave(time,signal):
    '''
    For fitting a sine wave to a signal trace

    Parameters
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
