import numpy as np
import numba
from numba.pycc import CC

cc = CC('model_func')

@cc.export("model_func", "float64[:](float64,float64[:],float64[:])")
def model_function(time, S, k):

    P = np.zeros(len(S))
    P[0] = -k[0]*S[0]*S[1]-k[1]*S[0]*S[3]
    P[1] = +k[1]*S[0]*S[3]-k[0]*S[0]*S[1]-k[2]*S[4]*S[1]-k[3]*S[4]*S[1]-k[4]*S[7]*S[1]
    P[2] = +k[0]*S[0]*S[1]
    P[3] = +k[1]*S[0]*S[3]-k[1]*S[0]*S[3]
    P[4] = +k[4]*S[7]*S[1]-k[2]*S[4]*S[1]-k[3]*S[4]*S[1]
    P[5] = +k[2]*S[4]*S[1]
    P[6] = +k[3]*S[4]*S[1]
    P[7] = +k[4]*S[7]*S[1]-k[4]*S[7]*S[1]

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
}

k = np.zeros(max(reactions.values())+1) # rate constants

S = np.zeros(len(species)) # initial concentrations

if __name__ == "__main__":
    cc.compile()
