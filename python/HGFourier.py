import numpy as np
import scipy.special

"""
    mu: the cosine of the output angle (elevation)
    mi: the cosine of the input angle (elevation)
    g: Henyey Greenstein scattering parameter [-1, 1]
    k: upper bound for fourier expansion
return
    coeff: the computed fourier coefficients (np.matrix)
"""


def HG(mu_i, mu_o, phi_d, g):
    sqrtPara = (1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o)
    cosTheta = mu_i * mu_o + np.cos(phi_d) * np.sqrt(sqrtPara)
    temp = (1.0 + g * g - 2.0 * g * cosTheta)
    return (1.0 - g *g) / (4.0 * np.pi * temp * np.sqrt(temp))


def HGFourier(mo, mi, g, k):
    # compute the A and B coefficients
    a = 1.0 + g * g - 2.0 * g * mi * mo
    b = -2.0 * g * np.sqrt((1.0 - mi * mi) * (1.0 - mo * mo))

    # The latter two integrals can be expressed in full elliptic integrals of the first, E(x), and second, K(x),
    absB = np.abs(b)
    arg = np.sqrt(2.0 * absB / (a + absB))
    K = scipy.special.ellipk(arg)
    E = scipy.special.ellipe(arg)
    sqrtAB = np.sqrt(a + absB)
    temp = (1.0 - g * g) * 0.5 / (np.pi * np.pi)

    # compute the first coefficient (CHECK PAPER (20))
    coeff0 = (E * temp * sqrtAB) / (a * a - b * b)

    coeff1 = 0.0

    if b != 0.0:
        # compute the seconde coefficient (CHECK PAPER (21))
        coeff1 = np.sign(b) * temp / (absB * sqrtAB) * (K - a / (a - absB) * E)

    coeff = np.zeros(k)

    m = max(k, 500)
    s = np.zeros(m + 1)

    # compute the ratio
    z = a / np.sqrt(a * a - b * b)
    delta = z / np.sqrt(z * z - 1.0)

    s[m] = (1.0 + 1.0 / (2.0 * m) - (1.0 + 3.0 * z) / (8.0 * m * m)) * np.sqrt((z - 1.0) / (z + 1.0))

    while True:
        m -= 1
        s[m] = (2.0 * m + 3.0) / (4.0 * (m + 1.0) * delta - (2.0 * m + 1.0) * s[m + 1])
        if m == 0:
            break

    C = 0.0
    if s[0] != 0.0:
        C = coeff1 / (coeff0 * s[0])

    coeff[0] = coeff0

    prod = coeff0 * C * 2.0

    for j in range(k - 1):
        prod *= s[j]
        if j % 2 == 0:
            coeff[j + 1] = prod
        else:
            coeff[j + 1] = prod * np.sign(g)

    return coeff





