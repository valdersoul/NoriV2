import numpy as np
import math
import scipy.special as sp
import matplotlib.pyplot as plt
from microfacet import microfacetNoExp, getAB, Bmax, microfacetFourierSeries, microfacetNoExpFourierSeries
from filon import filonIntegrate
from BSDFFourier import expcos_fseries, fseries_convolve, expcosCoefficientCount
from fourierBasic import fourierEval
import layerlab as ll


n = 12
phi = np.linspace(0.0,np.pi, n)

mu_o = -0.452930232232
mu_i = -0.024736727622
alpha = 0.1
etaC = complex(1.1, 0)

cLL = ll.microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, 10e-4)
cOur = microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, 10e-4)

print("Our coeff: " + str(len(cOur)))
print("ll coeff: " + str(len(cLL)))

plt.figure()
plt.plot(cOur, 'b')
plt.plot(cLL, 'r--')
plt.savefig("Coeff.pdf")
if len(cOur) == len(cLL):
    print("absolute error = " + str(np.sum(np.abs(np.array(cOur) - np.array(cLL)))))
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

