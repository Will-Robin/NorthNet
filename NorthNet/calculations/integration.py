def integrate_model(time, S0, k, model):
    '''
    For numerically integrating a model over time. Kind of old and not used in
    most things this module is used for.
    Should be modified and documented later.
    '''

    y = ode(execute_string_function).set_integrator('lsoda')
    y.set_initial_value(S0)
    y.set_f_params(k,model)

    calc = np.zeros([len(time),len(S0)])
    time_out = np.zeros([len(time)])

    d_time = np.hstack((0.,time))
    dt = [(y-x) for x,y in zip(d_time, d_time[1:])]

    for n in range(0,len(time)):
        calc[n] = y.integrate(y.t + dt[n])
        time_out[n] = y.t

    return time_out, calc

def integrate_model(function, S0, k, C, time, dt):

    y = ode(function).set_integrator('lsoda')
    y.set_initial_value(S0)
    y.set_f_params(k, C)

    calc = []
    time_arr = []
    while y.successful() and y.t <= time:
        calc.append(y.integrate(y.t+dt))
        time_arr.append(y.t+dt)

    calc = np.array(calc)
    time_arr = np.array(time_arr)

    return time_arr, calc
