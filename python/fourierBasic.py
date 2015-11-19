import numpy as np


def fourierEval(b, phi):
    y = np.zeros(len(phi))
    for i in range(len(phi)):
        curPhi = phi[i]
        for l in range(len(b)):
            y[i] += b[l] * np.cos(float(l) * curPhi)
    return y
