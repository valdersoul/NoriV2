import numpy as np
import math
import scipy.special as sp
import matplotlib.pyplot as plt
from microfacet import microfacetNoExp, getAB, Bmax
from filon import filonIntegrate
from BSDFFourier import expcos_fseries, fseries_convolve, expcosCoefficientCount
from fourierBasic import fourierEval
import layerlab as ll


n = 100
a = 0.0
b = np.pi
nCoeffb = 12
nEvals = 200
phi = np.linspace(a, b, n)
y  = np.zeros(n)
y2  = np.zeros(n)
mu_i = -0.6
mu_o = 0.5
alpha = 0.05
eta_i = 1.5
eta_o = 1.0
eta = eta_i / eta_o
etaC = complex(eta, 0)
relerr = 10e-4

for i in range(n):
    curPhi = phi[i]
    y[i] = microfacetNoExp(mu_o, mu_i, etaC, alpha, curPhi)


def f(a):
    return microfacetNoExp(mu_o, mu_i, etaC, alpha, a)
# compute the fourier approx
b = filonIntegrate(f, nCoeffb, nEvals, a, b)

A, B = getAB(mu_i, mu_o, etaC, alpha)

B_max = Bmax(n, relerr)
if B > B_max:
    A = A + B - B_max + math.log(sp.i0e(B) / sp.i0e(B_max))
    B = B_max

nCoeffa = expcosCoefficientCount(B, relerr)
print (nCoeffa)
a = expcos_fseries(A, B, nCoeffa)
c = fseries_convolve(a, nCoeffa, b, nCoeffb)
# c = ll.microfacetFourierSeries(mu_o, mu_i, etaC, alpha, nCoeff, 10e-4)

y = fourierEval(c, phi)

plt.figure()
plt.plot(phi, y, 'b')
plt.title("Used "+  str(nCoeffa) + " for the high frequency")
#plt.show()
plt.savefig("MicrofacetgraphAll.pdf")

