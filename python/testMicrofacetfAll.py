import numpy as np
import math
import scipy.special as sp
import matplotlib.pyplot as plt
from microfacet import microfacetNoExp, getAB, Bmax, microfacetFourierSeries, microfacetNoExpFourierSeries
from filon import filonIntegrate
from BSDFFourier import expcos_fseries, fseries_convolve, expcosCoefficientCount
from fourierBasic import fourierEval
import layerlab as ll


n = 200
phi = np.linspace(0.0, np.pi, n)

mu_o = -0.6
mu_i = 0.5
alpha = 0.05
etaC = complex(1.1, 0)


cOur = microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, 1e-3)


yOur = fourierEval(cOur, phi)

yLL = ll.microfacet(mu_o, mu_i, etaC, alpha, phi)

plt.figure()
plt.plot(phi, yLL, 'b')
x1,x2,y1,y2 = plt.axis()
plt.axis((y1,y2,0.0,1.0))
plt.axis((x1,x2,0.0,1.0))

plt.ylabel('Value')
plt.xlabel('Phi')
plt.title('Microfacet BSDF')
plt.savefig("MBSDF.png")

plt.figure()
plt.plot(phi, yOur, 'r')
x1,x2,y1,y2 = plt.axis()
plt.axis((y1,y2,0.0,1.0))
plt.axis((x1,x2,0.0,1.0))
plt.ylabel('Value')
plt.xlabel('Phi')
plt.title('Microfacet BSDF Fourier Approximation')
plt.savefig("MBSDFF.png")



