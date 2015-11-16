import numpy as np
import sys


# Search for the interior roots of P_{n+1}(x) using Newton's method.
def getRootOfLegendre(n, i):
    x = - np.cos((2.0 * float(i) + 1.0) / (2.0 * float(n) + 2.0) * np.pi)
    it = 0

    while True:
        it += 1
        if it > 20:
            print('did not converge after 20 iterations!')

        _, _, LNext, DNext = legendre(n + 1, x)
        step = LNext / DNext
        x -= step

        if np.abs(step) <= 4 * np.abs(x) * sys.float_info.epsilon:
            break
    return x


def legendre(l, x):
    """
        l: the level of the gaus lobatto function
        x:
    returns:
        nodes: the points on the x axis
        weights: the weights of the nodes
    """
    if l == 0:
        return 1, 0
    elif l == 1:
        return x, 1
    else:
        Lppred = 1
        Lpred = x
        Lcur = 0
        Dppred = 0
        Dpred = 1
        Dcur = 0
        for k in range(2, l):
            Lcur = ((2.0 * k - 1.0) * x * Lpred - (k - 1.0) * Lppred) / k
            Dcur = Dppred + (2.0 * k - 1.0) * Lpred
            Lppred = Lpred
            Lpred = Lcur
            Dppred = Dpred
            Dpred = Dcur

        Lnext = ((2.0 * l + 1.0) * x * Lpred - l * Lppred) / (l + 1.0)
        Dnext = Dppred + (2.0 * l + 1.0) * Lpred

        return Lcur, Dcur, Lnext - Lppred, Dnext - Dppred


def gausLobatto(n):
    """
        n: number of nodes
    returns:
        nodes: the points on the x axis
        weights: the weights of the nodes
    """
    nodes = np.zeros(n)
    weights = np.zeros(n)
    n -= 1
    if n < 2:
        print("n must be >= 1")

    nodes[0] = -1.0
    nodes[n] = 1.0

    weights[0] = 2.0 / (float(n) * (float(n) + 1.0))
    weights[n] = 2.0 / (float(n) * (float(n) + 1.0))

    m = int((n + 1) / 2)
    for i in range(1, m):
        x = getRootOfLegendre(n, i)

        Ln, _, _, _ = legendre(n + 1, x)
        weights[i] = (2.0 / ((n * (n + 1.0)) * Ln * Ln))
        weights[n - i] = (2.0 / ((n * (n + 1.0)) * Ln * Ln))
        nodes[i] = x
        nodes[n - i] = -x

    if (n % 2) == 0:
        L, _, _, _ = legendre(n + 1, 0.0)
        weights[int(n / 2)] = (2.0 / (L * L))
        nodes[int(n / 2)] = 0

    return nodes, weights
