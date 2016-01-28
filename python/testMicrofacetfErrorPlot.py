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
yLL = ll.microfacet(mu_o, mu_i, etaC, alpha, phi)

testRange = range(2, 150)
yError = []
for l in testRange:
    cOur = microfacetFourierSeries(mu_o, mu_i, etaC, alpha, l, 1e-3)
    yOur = fourierEval(cOur, phi)
    yError.append(np.sum(np.abs(yLL - yOur)))



plt.figure()
plt.plot(testRange, yError, 'b')
#x1,x2,y1,y2 = plt.axis()
#plt.axis((y1,y2,0.0,1.0))
#plt.axis((x1,x2,0.0,1.0))

plt.ylabel('Absolute error')
plt.xlabel('Number of coefficents')
plt.title('Error plot')
plt.savefig("MBSDFError.png")


