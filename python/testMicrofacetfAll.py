import numpy as np
import math
import scipy.special as sp
import matplotlib.pyplot as plt
from microfacet import microfacetNoExp, getAB, Bmax, microfacetFourierSeries, microfacetNoExpFourierSeries
from filon import filonIntegrate
from BSDFFourier import expcos_fseries, fseries_convolve, expcosCoefficientCount
from fourierBasic import fourierEval
import layerlab as ll


n = 1000
a = 0.0
b = np.pi
nCoeffb = 4
nEvals = 200
phi = np.linspace(a, b, n)
y  = np.zeros(n)
y2  = np.zeros(n)
mu_i = -0.6
mu_o = 0.5
alpha = 0.05
eta_i = 1.5c
eta_o = 1.0
eta = eta_i / eta_o
etaC = complex(eta, 0)
relerr = 10e-6

for i in range(n):
    curPhi = phi[i]
    y[i] = microfacetNoExp(mu_o, mu_i, etaC, alpha, curPhi)

quickTest = ll.microfacetFourierSeries(1.0, 1.0, etaC, alpha, n, 10e-4)
quickTest2 = microfacetFourierSeries(1.0, 1.0, etaC, alpha, n, 10e-4)
cOrig = ll.microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, 10e-4)
c = microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n)
cnoExp = microfacetNoExpFourierSeries(mu_o, mu_i, etaC, alpha, nCoeffb, np.pi)

print("Our coeff: " + str(len(c)))
print("ll coeff: " + str(len(cOrig)))

plt.figure()
plt.plot(cOrig, 'b')
plt.plot(c, 'r--')
plt.savefig("Coeff.pdf")

yOrig = fourierEval(cOrig, phi)
y = fourierEval(c, phi)

plt.figure()
plt.plot(phi, yOrig, 'b')
plt.plot(phi, y, 'r--')
plt.title("for the high frequency")
plt.savefig("MicrofacetgraphAll.pdf")

