import numpy as np
import scipy as sp
import math

def CM2S(M, normalizedFreq):
    R1 = 1.0
    RN = 1.0
    sqrt_R1_RN = math.sqrt(R1 * RN)
    N = M.shape[0] - 2

    U = np.eye(N + 2)
    U[0, 0] = 0.0
    U[-1, -1] = 0.0

    R = np.zeros(M.shape)
    R[0, 0] = R1
    R[-1, -1] = RN

    S21 = np.zeros(normalizedFreq.shape, dtype=np.complex)
    S11 = np.zeros(normalizedFreq.shape, dtype=np.complex)
    b = np.zeros((N + 2,), dtype=np.complex)
    b[0] = 1.0

    for k in np.arange(0, normalizedFreq.shape[0]):
        Z = normalizedFreq[k] * U - 1j * R + M

        Y = sp.linalg.lapack.zsysv(Z, b)[2]
        S11[k] = 1. + 2j * R1 * Y[0]
        S21[k] = -2j * sqrt_R1_RN * Y[-1]

    return S11, S21