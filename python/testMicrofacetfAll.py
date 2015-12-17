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
eta_i = 1.5
eta_o = 1.0
eta = eta_i / eta_o
etaC = complex(eta, 0)
relerr = 10e-6

for i in range(n):
    curPhi = phi[i]
    y[i] = microfacetNoExp(mu_o, mu_i, etaC, alpha, curPhi)


cLL = ll.microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, 10e-4)
cOur = microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, 10e-4)

print("Our coeff: " + str(len(cOur)))
print("ll coeff: " + str(len(cLL)))

plt.figure()
plt.plot(cOur, 'b')
plt.plot(cLL, 'r--')
plt.savefig("Coeff.pdf")
if len(cOur) == len(cLL):
    print("absolute error = " + str(np.sum(np.abs(cOur - cLL ))))
else:
    print("len(expcos_coeffs) = " + str(len(cOur)))
    print("len(llexpcos_coeffs) = " + str(len(cLL )))


yOur = fourierEval(cOur, phi)

yLL = ll.microfacet(mu_o, mu_i, etaC, alpha, phi)

plt.figure()
plt.plot(phi, yOur, 'b')
plt.plot(phi, yLL, 'r--')
plt.ylabel('Value')
plt.xlabel('Phi_d')
plt.title('Microfacet fourier approximation')
plt.legend(['Fourier approximation', 'Orignal function'])
plt.savefig("MicrofacetgraphAll.png")

