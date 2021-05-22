from scipy.integrate import ode
import numpy as np

def integrate_model(function, S0, k, time, dt):

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
