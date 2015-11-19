import numpy as np6
import matplotlib.pyplot as plt
from microfacet import microfacetNoExp
from filon import filonIntegrate


n = 100
a = 0.0
b = np.pi
nCoeff = 6
nEvals = 200
phi = np.linspace(a, b, n)
y  = np.zeros(n)
y2  = np.zeros(n)
mu_i = -0.6
mu_o = 0.5
alpha = 0.05
eta = 1.5

for i in range(n):
    curPhi = phi[i]
    etaC = complex(eta, 0)
    y[i] = microfacetNoExp(mu_o, mu_i, etaC, alpha, curPhi)


def f(a):
    return microfacetNoExp(mu_o, mu_i, etaC, alpha, a)
# compute the fourier approx
b = filonIntegrate(f, nCoeff, nEvals, a, b)

for i in range(n):
    curPhi = phi[i]
    for l in range(nCoeff):
        y2[i] += b[l] * np.cos(float(l) * curPhi)

plt.figure()
plt.plot(phi, y, 'b')
plt.plot(phi, y2, 'r--')
err = np.sum(np.abs(y - y2)) / n
plt.title("Error " + str(err) + " with " + str(nCoeff) + " coefficients")
#plt.show()
plt.savefig("Microfacetgraph.pdf")

