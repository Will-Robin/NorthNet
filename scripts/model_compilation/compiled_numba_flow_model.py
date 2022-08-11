import numpy as np
import numba
from numba.pycc import CC

cc = CC('model_func')

@cc.export("model_func", "float64[:](float64,float64[:],float64[:])")
def model_function(time, S, k):

    P = np.zeros(len(S))

    F_in = np.array(
    [[0.006291058,0.006291058,0.006291058,0.006291058,0.006291058,0.000884192,
      0.000884192,0.000884192,0.000884192],
     [0.000347222,0.000347222,0.000347222,0.000347222,0.000347222,0.000347222,
      0.000347222,0.000347222,0.000347222],
     [0.001139497,0.001139497,0.001139497,0.001139497,0.001139497,0.006546364,
      0.006546364,0.006546364,0.006546364],
     [0.000000000,0.000000000,0.000000000,0.000000000,0.000000000,0.000000000,
      0.000000000,0.000000000,0.000000000]]    )

    flow_time = np.array(
    [0.000000000,200.000000000,400.000000000,600.000000000,800.000000000,
     1000.000000000,1200.000000000,1400.000000000,1600.000000000]    )

    total_flow = np.array(
    [[0.008125000,0.008125000,0.008125000,0.008125000,0.008125000,0.008125000,
      0.008125000,0.008125000,0.008125000]]    )

    i = np.abs(flow_time - time).argmin()

    P[0] = +(F_in[0,i]*2.0)-k[0]*S[0]*S[1]-k[1]*S[0]*S[3]-S[0]*total_flow[0,i]
    P[1] = +k[1]*S[0]*S[3]-k[0]*S[0]*S[1]-k[2]*S[4]*S[1]-k[3]*S[4]*S[1]-k[4]*S[7]*S[1]-S[1]*total_flow[0,i]
    P[2] = +k[0]*S[0]*S[1]-S[2]*total_flow[0,i]
    P[3] = +k[1]*S[0]*S[3]+(F_in[1,i]*0.12)-k[1]*S[0]*S[3]-S[3]*total_flow[0,i]
    P[4] = +k[4]*S[7]*S[1]-k[2]*S[4]*S[1]-k[3]*S[4]*S[1]-S[4]*total_flow[0,i]
    P[5] = +k[2]*S[4]*S[1]-S[5]*total_flow[0,i]
    P[6] = +k[3]*S[4]*S[1]-S[6]*total_flow[0,i]
    P[7] = +k[4]*S[7]*S[1]+(F_in[2,i]*55.5)-k[4]*S[7]*S[1]-S[7]*total_flow[0,i]

    return P

def wrapper_function(time, S, k):
    return model_function(time, S, k)

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

inputs = {
'O=C(CO)CO':2.0,
'[OH-]':0.12,
'O':55.5,
}

k = np.zeros(max(reactions.values())+1) # rate constants

S = np.zeros(len(species)) # initial concentrations

if __name__ == "__main__":
    cc.compile()
