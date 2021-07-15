import numpy as np
from scipy.integrate import ode

def integrate_model(function, S0, k, time, dt):
    '''
    Wrapper for scipy ode
    uses old SciPy API (as of June 2021)

    function: function for integration
        function(time, variables, rate constants)

    S0: numpy array
        Initial state of the reaction
    k: numpy array
        Rate constants
    time: float
        The simulation will run from [0,time)
    dt: float
        time step

    time_arr: numpy array
        Time values from the integration

    calc: numpy array
        Calculated time progresses
    '''


    y = ode(function).set_integrator('lsoda')
    y.set_initial_value(S0)
    y.set_f_params(k)

    calc = []
    time_arr = []
    while y.successful() and y.t <= time:
        calc.append(y.integrate(y.t+dt))
        time_arr.append(y.t+dt)

    calc = np.array(calc)
    time_arr = np.array(time_arr)

    return time_arr, calc
