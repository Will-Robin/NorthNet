import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

from NorthNet.simulation.integration import integrate_model

import model_func as currentModel

def run_integration(S0, rate_constants, max_time, time_step):

    integrator = ode(currentModel.model_func).set_integrator("lsoda")
    integrator.set_initial_value(S0)
    integrator.set_f_params(rate_constants)

    calc = []
    time_arr = []
    while integrator.successful() and integrator.t <= max_time:
        calc.append(integrator.integrate(integrator.t + time_step))
        time_arr.append(integrator.t + time_step)

    calc = np.array(calc)
    time = np.array(time_arr)

    compound_traces = calc.T

    return time, compound_traces


def main(S, rate_constants, max_time, time_step):

    time, compound_traces = run_integration(S, rate_constants, max_time, time_step)

    fig, ax = plt.subplots()

    for x in range(0, len(compound_traces)):
        ax.plot(time, compound_traces[x])

    ax.set_xlabel("time/ s")
    ax.set_ylabel("conc./ M")

    plt.show()
    plt.close()


if __name__ == "__main__":

    species = {
    'O=C(CO)CO':0,
    'OC=C(O)CO':1,
    'O=C[C@@](O)(CO)C(O)(CO)CO':2,
    '[OH-]':3,
    'O=C[C@H](O)CO':4,
    'O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO':5,
    'O=C[C@@](O)(CO)[C@@H](O)[C@H](O)CO':6,
    'O':7,
    }

    reactions = {
    'O=C(CO)CO.OC=C(O)CO>>O=C[C@@](O)(CO)C(O)(CO)CO':0,
    'O=C(CO)CO.[OH-]>>OC=C(O)CO.[OH-]':1,
    'O=C[C@H](O)CO.OC=C(O)CO>>O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO':2,
    'O=C[C@H](O)CO.OC=C(O)CO>>O=C[C@@](O)(CO)[C@@H](O)[C@H](O)CO':3,
    'O.OC=C(O)CO>>O.O=C[C@H](O)CO':4,
    }

    k = np.zeros(max(reactions.values())+1) # rate constants

    S = np.zeros(len(species)) # initial concentrations

    S0 = np.zeros(S.shape)

    S0[0] = 2.0
    S0[7] = 0.12
    S0[3] = 0.12

    rate_constants = np.ones(k.shape) + 1000

    time_offset = 4318.0
    lead_in_time = 0.0

    main(S0, rate_constants, 10000.0, 1.0)
