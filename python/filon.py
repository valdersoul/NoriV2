import numpy as np


def filon(phi, f, lmax, output):

    h = phi[1] - phi[0]
    invH = 1 / h
    output[0] += (1.0 / (6.0 * np.pi)) * h * (f[0] + 4 * f[1] + f[2])

    cosPhi0Prev = np.cos(phi[0])
    cosPhi0Cur = 1.0
    cosPhi1Prev = np.cos(phi[1])
    cosPhi1Cur = 1.0
    sinPhi0Prev = -np.sin(phi[0])
    sinPhi0Cur = 0.0
    sinPhi1Prev = -np.sin(phi[1])
    sinPhi1Cur = 0.0
    twoCosPhi0 = 2.0 * cosPhi0Prev
    twoCosPhi1 = 2.0 * cosPhi1Prev

    term0 = 3 * f[0] - 4 * f[1] + f[2]
    term1 = f[0] - 4 * f[1] + 3 * f[2]
    term2 = 4 * (f[0] - 2 * f[1] + f[2])

    for l in range(1, lmax):
        cosPhi0Next = twoCosPhi0 * cosPhi0Cur - cosPhi0Prev
        cosPhi1Next = twoCosPhi1 * cosPhi1Cur - cosPhi1Prev
        sinPhi0Next = twoCosPhi0 * sinPhi0Cur - sinPhi0Prev
        sinPhi1Next = twoCosPhi1 * sinPhi1Cur - sinPhi1Prev

        invL = 1.0 / float(l)
        invL2H = invH * invL * invL
        invL3H2 = invL2H * invL * invH

        output[l] += (2 * (1.0 / np.pi)) * ((invL2H * (term0 * cosPhi0Next + term1 * cosPhi1Next) +
                                             invL3H2 * term2 * (sinPhi0Next - sinPhi1Next) +
                                             invL * (f[2] * sinPhi1Next - f[0] * sinPhi0Next)))

        cosPhi0Prev = cosPhi0Cur
        cosPhi0Cur = cosPhi0Next
        cosPhi1Prev = cosPhi1Cur
        cosPhi1Cur = cosPhi1Next
        sinPhi0Prev = sinPhi0Cur
        sinPhi0Cur = sinPhi0Next
        sinPhi1Prev = sinPhi1Cur
        sinPhi1Cur = sinPhi1Next

    return output


def filonIntegrate(f, nCoeffs, nEvals, a, b):
    coeffs = np.zeros(nCoeffs)
    if nEvals % 2 == 0:
        nEvals += 1

    value = np.zeros(3)
    phi = np.zeros(2)
    delta = (b - a) / (nEvals - 1)
    phi[0] = a
    value[0] = f(a)
    for i in range(int((nEvals - 1) / 2)):
        phi[1] = phi[0] + 2 * delta
        value[1] = f(phi[0] + delta)
        value[2] = f(phi[1])

        coeffs = filon(phi, value, nCoeffs, coeffs)

        value[0] = value[2]
        phi[0] = phi[1]

    return coeffs
