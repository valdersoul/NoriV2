import numpy as np
from gaussLobatto import gausLobatto

def integrand(x):
    return 6.3 * (x**2 - x**6 - x**8)

print (np.sum(integrand(np.linspace(-1, 1, 100)) * 1.0/50.0))
n = 6
mu, w = gausLobatto(n)

print (np.sum(integrand(mu) * w))
