import numpy as np
import sys
import scipy.special as sp


def fseries_convolve(a, ka, b, kb, c):
    for i in range(ka + kb - 1):
        sum = 0.0
        for j in range(np.min(kb, ka - i)):
            sum += b[j] * a[i + j]

        for j in range(np.max(0, i - ka + 1), np.min(kb, i + ka)):
            sum += b[j] * a[abs(i - j)]
        if i < kb:
            sum += b[i] * a[0]
        if i == 0:
            sum = 0.5 * (sum + a[0] * b[0])
        c[i] = 0.5 * sum
    return c


def mbessel_ratio(k, B):
    eps = sys.float_info.epsilon
    invTwoB = 2.0 / B
    i = k
    D = 1.0 / (invTwoB * i)
    i += 1
    Cd = D
    C = Cd

    while np.abs(Cd) > eps * np.abs(C):
        coeff = invTwoB * i
        i += 1
        D = 1 / (D + coeff)
        Cd *= coeff * D - 1.0
        C += Cd

    return C


def expcos_fseries(A, B, m):
    coeffs = np.zeros(m + 1)
    # /* Determine the last ratio and work downwards */
    coeffs[m] = mbessel_ratio(B, m)
    for i in range(m - 1, 1):
        coeffs[i] = B / (2 * i + B * coeffs[i + 1])
        # /* Evaluate the exponentially scaled I0 and correct scaling */
    coeffs[0] = np.exp(A + B) * sp.i0e(B)
    # /* Apply the ratios & factor of two upwards */
    prod = 2 * coeffs[0]
    for i in range(m + 1):
        prod *= coeffs[i]
        coeffs[i] = prod

    return coeffs

