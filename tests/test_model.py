import numpy as np
import numba


@numba.jit(
    numba.float64[:](numba.float64, numba.float64[:], numba.float64[:]),
    locals={"P": numba.float64[:], "F": numba.float64[:, :], "I": numba.float64[:]},
    nopython=True,
)
def model_function(time, S, k):

    P = np.zeros(len(S))

    F = np.array(
        [
            [0.000000000],
            [0.000000000],
            [0.000000000],
            [0.000000000],
            [0.000000000],
            [0.000000000],
            [0.008333333],
        ]
    )

    idx = np.abs(F[0] - time).argmin()

    I = F[1:-1, idx]

    sigma_flow = F[-1, idx]

    P[0] = -k[0] * S[0] * S[1] - k[16] * S[0] * S[2]
    P[1] = (
        +k[0] * S[0] * S[1]
        + k[4] * S[8] * S[1]
        + k[9] * S[6] * S[1]
        + k[13] * S[8] * S[1]
        + k[20] * S[13] * S[1]
        + k[22] * S[25] * S[1]
        + k[25] * S[25] * S[1]
        + k[27] * S[4] * S[1]
        + k[28] * S[12] * S[1]
        + k[32] * S[11] * S[1]
        + k[36] * S[7] * S[1]
        + k[38] * S[11] * S[1]
        + k[43] * S[33] * S[1]
        + k[50] * S[19] * S[1]
        + k[54] * S[27] * S[1]
        - k[0] * S[0] * S[1]
        - k[4] * S[8] * S[1]
        - k[9] * S[6] * S[1]
        - k[13] * S[8] * S[1]
        - k[18] * S[10] * S[22] * S[1]
        - k[20] * S[13] * S[1]
        - k[22] * S[25] * S[1]
        - k[25] * S[25] * S[1]
        - k[27] * S[4] * S[1]
        - k[28] * S[12] * S[1]
        - k[32] * S[11] * S[1]
        - k[36] * S[7] * S[1]
        - k[38] * S[11] * S[1]
        - k[43] * S[33] * S[1]
        - k[50] * S[19] * S[1]
        - k[54] * S[27] * S[1]
        - k[55] * S[10] * S[12] * S[1]
    )
    P[2] = (
        +k[0] * S[0] * S[1]
        + k[6] * S[8]
        + k[27] * S[4] * S[1]
        + k[43] * S[33] * S[1]
        + k[47] * S[17]
        + k[49] * S[16]
        + k[52] * S[20]
        + k[53] * S[25]
        - k[1] * S[2] * S[3]
        - k[3] * S[7] * S[2]
        - k[8] * S[2] * S[10]
        - k[11] * S[4] * S[2]
        - k[12] * S[7] * S[2]
        - k[15] * S[4] * S[2]
        - k[16] * S[0] * S[2]
        - k[19] * S[4] * S[2]
        - k[56] * S[2] * S[10]
    )
    P[3] = (
        +k[1] * S[2]
        + k[14] * S[18]
        + k[23] * S[28]
        + k[26] * S[30]
        + k[31] * S[26]
        + k[33] * S[32]
        + k[39] * S[30]
        + k[40] * S[18]
        + k[44] * S[26]
        + k[51] * S[18]
        - k[1] * S[2] * S[3]
        - k[14] * S[18] * S[3]
        - k[23] * S[28] * S[3]
        - k[26] * S[30] * S[3]
        - k[31] * S[26] * S[3]
        - k[33] * S[32] * S[3]
        - k[39] * S[30] * S[3]
        - k[40] * S[18] * S[3]
        - k[44] * S[26] * S[3]
        - k[51] * S[18] * S[3]
    )
    P[4] = (
        +k[1] * S[2]
        + k[24] * S[6]
        + k[37] * S[5] * S[10]
        + k[49] * S[16]
        + k[52] * S[20]
        + k[53] * S[25]
        - k[2] * S[4] * S[5]
        - k[11] * S[4] * S[2]
        - k[15] * S[4] * S[2]
        - k[19] * S[4] * S[2]
        - k[27] * S[4] * S[1]
    )
    P[5] = (
        +k[24] * S[6]
        + k[36] * S[7] * S[1]
        + k[41] * S[12]
        + k[45] * S[22]
        - k[2] * S[4] * S[5]
        - k[7] * S[7] * S[5]
        - k[17] * S[7] * S[5]
        - k[37] * S[5] * S[10]
        - k[42] * S[5] * S[10]
    )
    P[6] = +k[2] * S[4] * S[5] - k[9] * S[6] * S[1] - k[24] * S[6]
    P[7] = (
        +k[6] * S[8]
        + k[30] * S[11]
        + k[41] * S[12]
        + k[45] * S[22]
        + k[47] * S[17]
        + k[48] * S[15]
        - k[3] * S[7] * S[2]
        - k[7] * S[7] * S[5]
        - k[12] * S[7] * S[2]
        - k[17] * S[7] * S[5]
        - k[36] * S[7] * S[1]
    )
    P[8] = (
        +k[3] * S[7] * S[2]
        + k[51] * S[18]
        - k[4] * S[8] * S[1]
        - k[6] * S[8]
        - k[13] * S[8] * S[1]
    )
    P[9] = +k[4] * S[8] * S[1] - k[5] * S[9] * S[10] - k[34] * S[9] * S[10]
    P[10] = (
        -k[5] * S[9] * S[10]
        - k[8] * S[2] * S[10]
        - k[10] * S[14] * S[10]
        - k[18] * S[10] * S[22] * S[1]
        - k[21] * S[26] * S[10]
        - k[29] * S[26] * S[10]
        - k[34] * S[9] * S[10]
        - k[35] * S[31] * S[10]
        - k[37] * S[5] * S[10]
        - k[42] * S[5] * S[10]
        - k[46] * S[26] * S[10]
        - k[55] * S[10] * S[12] * S[1]
        - k[56] * S[2] * S[10]
        - k[57] * S[18] * S[10]
    )
    P[11] = (
        +k[5] * S[9] * S[10]
        + k[26] * S[30]
        - k[30] * S[11]
        - k[32] * S[11] * S[1]
        - k[38] * S[11] * S[1]
    )
    P[12] = (
        +k[7] * S[7] * S[5]
        - k[28] * S[12] * S[1]
        - k[41] * S[12]
        - k[55] * S[10] * S[12] * S[1]
    )
    P[13] = +k[8] * S[2] * S[10] + k[31] * S[26] - k[20] * S[13] * S[1]
    P[14] = +k[9] * S[6] * S[1] - k[10] * S[14] * S[10]
    P[15] = +k[10] * S[14] * S[10] + k[33] * S[32] - k[48] * S[15]
    P[16] = +k[11] * S[4] * S[2] + k[34] * S[9] * S[10] - k[49] * S[16]
    P[17] = +k[12] * S[7] * S[2] + k[46] * S[26] * S[10] - k[47] * S[17]
    P[18] = (
        +k[13] * S[8] * S[1]
        + k[50] * S[19] * S[1]
        + k[54] * S[27] * S[1]
        - k[14] * S[18] * S[3]
        - k[40] * S[18] * S[3]
        - k[51] * S[18] * S[3]
        - k[57] * S[18] * S[10]
    )
    P[19] = (
        +k[14] * S[18]
        + k[29] * S[26] * S[10]
        + k[35] * S[31] * S[10]
        - k[50] * S[19] * S[1]
    )
    P[20] = +k[15] * S[4] * S[2] - k[52] * S[20]
    P[21] = +k[16] * S[0] * S[2]
    P[22] = (
        +k[17] * S[7] * S[5]
        + k[44] * S[26]
        - k[18] * S[10] * S[22] * S[1]
        - k[45] * S[22]
    )
    P[23] = +k[18] * S[10] * S[22] * S[1]
    P[24] = +k[18] * S[10] * S[22] * S[1] + k[55] * S[10] * S[12] * S[1]
    P[25] = (
        +k[19] * S[4] * S[2]
        + k[39] * S[30]
        - k[22] * S[25] * S[1]
        - k[25] * S[25] * S[1]
        - k[53] * S[25]
    )
    P[26] = (
        +k[20] * S[13] * S[1]
        + k[28] * S[12] * S[1]
        - k[21] * S[26] * S[10]
        - k[29] * S[26] * S[10]
        - k[31] * S[26] * S[3]
        - k[44] * S[26] * S[3]
        - k[46] * S[26] * S[10]
    )
    P[27] = +k[21] * S[26] * S[10] + k[40] * S[18] - k[54] * S[27] * S[1]
    P[28] = +k[22] * S[25] * S[1] - k[23] * S[28] * S[3]
    P[29] = +k[23] * S[28]
    P[30] = (
        +k[25] * S[25] * S[1]
        + k[38] * S[11] * S[1]
        - k[26] * S[30] * S[3]
        - k[39] * S[30] * S[3]
    )
    P[31] = +k[30] * S[11] + k[48] * S[15] - k[35] * S[31] * S[10]
    P[32] = +k[32] * S[11] * S[1] - k[33] * S[32] * S[3]
    P[33] = +k[42] * S[5] * S[10] - k[43] * S[33] * S[1]
    P[34] = +k[55] * S[10] * S[12] * S[1]
    P[35] = +k[56] * S[2] * S[10]
    P[36] = +k[57] * S[18] * S[10]
    return P


def wrapper_function(time, S, k):
    return model_function(time, S, k)


species = {
    "O=C(CO)CO": 0,
    "[OH-]": 1,
    "OC=C(O)CO": 2,
    "O": 3,
    "O=C[C@H](O)CO": 4,
    "OC=CO": 5,
    "O=C[C@@H](O)[C@@H](O)[C@H](O)CO": 6,
    "O=CCO": 7,
    "O=C(CO)[C@H](O)[C@H](O)CO": 8,
    "OC=C(O)[C@H](O)[C@H](O)CO": 9,
    "C=O": 10,
    "O=C([C@@H](O)CO)[C@H](O)[C@H](O)CO": 11,
    "O=C[C@@H](O)[C@H](O)CO": 12,
    "O=C(CO)[C@H](O)CO": 13,
    "OC=C(O)[C@@H](O)[C@H](O)CO": 14,
    "O=C([C@@H](O)CO)[C@@H](O)[C@H](O)CO": 15,
    "O=C[C@@](O)(CO)[C@H](O)[C@H](O)CO": 16,
    "O=C[C@@](O)(CO)[C@H](O)CO": 17,
    "OCC(O)=C(O)[C@H](O)CO": 18,
    "O=C([C@@H](O)CO)[C@H](O)CO": 19,
    "O=C[C@@](O)(CO)[C@@H](O)[C@H](O)CO": 20,
    "O=C[C@@](O)(CO)C(O)(CO)CO": 21,
    "O=C[C@H](O)[C@H](O)CO": 22,
    "OC[C@H](O)[C@H](O)CO": 23,
    "O=CO": 24,
    "O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO": 25,
    "OC=C(O)[C@H](O)CO": 26,
    "O=C([C@H](O)CO)[C@H](O)CO": 27,
    "OC=C(O)[C@H](O)[C@H](O)[C@H](O)CO": 28,
    "O=C[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO": 29,
    "OCC(O)=C(O)[C@H](O)[C@H](O)CO": 30,
    "OC=C(O)[C@@H](O)CO": 31,
    "OC[C@H](O)C(O)=C(O)[C@H](O)CO": 32,
    "O=C[C@@H](O)CO": 33,
    "OC[C@@H](O)[C@H](O)CO": 34,
    "O=CC(O)(CO)CO": 35,
    "O=C(CO)[C@](O)(CO)[C@H](O)CO": 36,
}

reactions = {
    "O=C(CO)CO.[OH-]>>OC=C(O)CO.[OH-]": 0,
    "O.OC=C(O)CO>>O.O=C[C@H](O)CO": 1,
    "O=C[C@H](O)CO.OC=CO>>O=C[C@@H](O)[C@@H](O)[C@H](O)CO": 2,
    "O=CCO.OC=C(O)CO>>O=C(CO)[C@H](O)[C@H](O)CO": 3,
    "O=C(CO)[C@H](O)[C@H](O)CO.[OH-]>>OC=C(O)[C@H](O)[C@H](O)CO.[OH-]": 4,
    "C=O.OC=C(O)[C@H](O)[C@H](O)CO>>O=C([C@@H](O)CO)[C@H](O)[C@H](O)CO": 5,
    "O=C(CO)[C@H](O)[C@H](O)CO>>O=CCO.OC=C(O)CO": 6,
    "O=CCO.OC=CO>>O=C[C@@H](O)[C@H](O)CO": 7,
    "C=O.OC=C(O)CO>>O=C(CO)[C@H](O)CO": 8,
    "O=C[C@@H](O)[C@@H](O)[C@H](O)CO.[OH-]>>OC=C(O)[C@@H](O)[C@H](O)CO.[OH-]": 9,
    "C=O.OC=C(O)[C@@H](O)[C@H](O)CO>>O=C([C@@H](O)CO)[C@@H](O)[C@H](O)CO": 10,
    "O=C[C@H](O)CO.OC=C(O)CO>>O=C[C@@](O)(CO)[C@H](O)[C@H](O)CO": 11,
    "O=CCO.OC=C(O)CO>>O=C[C@@](O)(CO)[C@H](O)CO": 12,
    "O=C(CO)[C@H](O)[C@H](O)CO.[OH-]>>OCC(O)=C(O)[C@H](O)CO.[OH-]": 13,
    "O.OCC(O)=C(O)[C@H](O)CO>>O.O=C([C@@H](O)CO)[C@H](O)CO": 14,
    "O=C[C@H](O)CO.OC=C(O)CO>>O=C[C@@](O)(CO)[C@@H](O)[C@H](O)CO": 15,
    "O=C(CO)CO.OC=C(O)CO>>O=C[C@@](O)(CO)C(O)(CO)CO": 16,
    "O=CCO.OC=CO>>O=C[C@H](O)[C@H](O)CO": 17,
    "C=O.O=C[C@H](O)[C@H](O)CO.[OH-]>>O=CO.OC[C@H](O)[C@H](O)CO": 18,
    "O=C[C@H](O)CO.OC=C(O)CO>>O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO": 19,
    "O=C(CO)[C@H](O)CO.[OH-]>>OC=C(O)[C@H](O)CO.[OH-]": 20,
    "C=O.OC=C(O)[C@H](O)CO>>O=C([C@H](O)CO)[C@H](O)CO": 21,
    "O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO.[OH-]>>OC=C(O)[C@H](O)[C@H](O)[C@H](O)CO.[OH-]": 22,
    "O.OC=C(O)[C@H](O)[C@H](O)[C@H](O)CO>>O.O=C[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO": 23,
    "O=C[C@@H](O)[C@@H](O)[C@H](O)CO>>O=C[C@H](O)CO.OC=CO": 24,
    "O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO.[OH-]>>OCC(O)=C(O)[C@H](O)[C@H](O)CO.[OH-]": 25,
    "O.OCC(O)=C(O)[C@H](O)[C@H](O)CO>>O.O=C([C@@H](O)CO)[C@H](O)[C@H](O)CO": 26,
    "O=C[C@H](O)CO.[OH-]>>OC=C(O)CO.[OH-]": 27,
    "O=C[C@@H](O)[C@H](O)CO.[OH-]>>OC=C(O)[C@H](O)CO.[OH-]": 28,
    "C=O.OC=C(O)[C@H](O)CO>>O=C([C@@H](O)CO)[C@H](O)CO": 29,
    "O=C([C@@H](O)CO)[C@H](O)[C@H](O)CO>>O=CCO.OC=C(O)[C@@H](O)CO": 30,
    "O.OC=C(O)[C@H](O)CO>>O.O=C(CO)[C@H](O)CO": 31,
    "O=C([C@@H](O)CO)[C@H](O)[C@H](O)CO.[OH-]>>OC[C@H](O)C(O)=C(O)[C@H](O)CO.[OH-]": 32,
    "O.OC[C@H](O)C(O)=C(O)[C@H](O)CO>>O.O=C([C@@H](O)CO)[C@@H](O)[C@H](O)CO": 33,
    "C=O.OC=C(O)[C@H](O)[C@H](O)CO>>O=C[C@@](O)(CO)[C@H](O)[C@H](O)CO": 34,
    "C=O.OC=C(O)[C@@H](O)CO>>O=C([C@@H](O)CO)[C@H](O)CO": 35,
    "O=CCO.[OH-]>>OC=CO.[OH-]": 36,
    "C=O.OC=CO>>O=C[C@H](O)CO": 37,
    "O=C([C@@H](O)CO)[C@H](O)[C@H](O)CO.[OH-]>>OCC(O)=C(O)[C@H](O)[C@H](O)CO.[OH-]": 38,
    "O.OCC(O)=C(O)[C@H](O)[C@H](O)CO>>O.O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO": 39,
    "O.OCC(O)=C(O)[C@H](O)CO>>O.O=C([C@H](O)CO)[C@H](O)CO": 40,
    "O=C[C@@H](O)[C@H](O)CO>>O=CCO.OC=CO": 41,
    "C=O.OC=CO>>O=C[C@@H](O)CO": 42,
    "O=C[C@@H](O)CO.[OH-]>>OC=C(O)CO.[OH-]": 43,
    "O.OC=C(O)[C@H](O)CO>>O.O=C[C@H](O)[C@H](O)CO": 44,
    "O=C[C@H](O)[C@H](O)CO>>O=CCO.OC=CO": 45,
    "C=O.OC=C(O)[C@H](O)CO>>O=C[C@@](O)(CO)[C@H](O)CO": 46,
    "O=C[C@@](O)(CO)[C@H](O)CO>>O=CCO.OC=C(O)CO": 47,
    "O=C([C@@H](O)CO)[C@@H](O)[C@H](O)CO>>O=CCO.OC=C(O)[C@@H](O)CO": 48,
    "O=C[C@@](O)(CO)[C@H](O)[C@H](O)CO>>O=C[C@H](O)CO.OC=C(O)CO": 49,
    "O=C([C@@H](O)CO)[C@H](O)CO.[OH-]>>OCC(O)=C(O)[C@H](O)CO.[OH-]": 50,
    "O.OCC(O)=C(O)[C@H](O)CO>>O.O=C(CO)[C@H](O)[C@H](O)CO": 51,
    "O=C[C@@](O)(CO)[C@@H](O)[C@H](O)CO>>O=C[C@H](O)CO.OC=C(O)CO": 52,
    "O=C(CO)[C@H](O)[C@H](O)[C@H](O)CO>>O=C[C@H](O)CO.OC=C(O)CO": 53,
    "O=C([C@H](O)CO)[C@H](O)CO.[OH-]>>OCC(O)=C(O)[C@H](O)CO.[OH-]": 54,
    "C=O.O=C[C@@H](O)[C@H](O)CO.[OH-]>>O=CO.OC[C@@H](O)[C@H](O)CO": 55,
    "C=O.OC=C(O)CO>>O=CC(O)(CO)CO": 56,
    "C=O.OCC(O)=C(O)[C@H](O)CO>>O=C(CO)[C@](O)(CO)[C@H](O)CO": 57,
}

inputs = {}

k = np.zeros(max(reactions.values()) + 1)  # rate constants

S = np.zeros(len(species))  # initial concentrations

C = np.zeros(len(inputs))  # input concentrations

time_offset = 9000.0
lead_in_time = 8130.25
