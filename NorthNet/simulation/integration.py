import numpy as np
from scipy.integrate import ode

def integrate_model(
                    function,
                    initial_species_concentrations,
                    k,
                    time,
                    time_interval
                    ):
    '''
    Wrapper for scipy ode
    uses old SciPy API (as of June 2021)

    Parameters
    ----------
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

    Returns
    -------
    time_arr: numpy array
        Time values from the integration

    calc: numpy array
        Calculated time progresses
    '''

    integrator = ode(function).set_integrator('lsoda')
    integrator.set_initial_value(initial_species_concentrations)
    integrator.set_f_params(k)

    calc = []
    time_arr = []
    while integrator.successful() and integrator.t <= time:
        calc.append(integrator.integrate(integrator.t+time_interval))
        time_arr.append(integrator.t+time_interval)

    calc = np.array(calc)
    time_arr = np.array(time_arr)

    return time_arr, calc

