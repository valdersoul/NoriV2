from HGFourier import *
import matplotlib.pyplot as plt
import numpy as np
import layerlab as ll


def angleSphere(mu, mu2, phi, phi2):
    v1 = np.matrix([np.sin(mu) * np.cos(phi),
                    np.sin(mu) * np.sin(phi),
                    np.cos(phi)
                    ])
    v2 = np.matrix([np.sin(mu2) * np.cos(phi2),
                    np.sin(mu2) * np.sin(phi2),
                    np.cos(phi2)
                    ]).T
    return v1.dot(v2)[0, 0]


def HGAPProx(mui, muo, phii, phio, g, k):
    result = 0.0
    coeff = HGFourier(np.cos(mui), np.cos(muo), g, k)
    relErr = 0.001
    #coeff = ll.hgFourierSeries(np.cos(muo), np.cos(mui), g, k, relErr)
    for l in range(len(coeff)):
        result += coeff[l] * np.cos(l * (phii - phio))
    return result


mu_o = np.pi / 4.0
mu_i = np.pi / 3.0
phi_o = 0.0
phi_i = np.pi

g = 0.75
cosGamma = angleSphere(np.cos(mu_o), np.cos(mu_i), abs(np.cos(phi_o) - np.cos(phi_i)), g)

realValue = HG(np.cos(mu_i), np.cos(mu_o), phi_i - phi_o, g)
#realValue = ll.hg(mu_o=[np.cos(mu_o)], mu_i=[np.cos(mu_i)], g=[g], phi_d=[(phi_i) - (phi_o)])
errList = []
xRange = range(3, 50)
for k in xRange:
    approxValue = HGAPProx(mu_o, mu_i, phi_o, phi_i, g, k)
    errList.append(abs(realValue - approxValue))

plt.figure()
plt.plot(xRange, errList)
plt.ylabel('abs(HG - HGAPProx_k)')
plt.xlabel('k - #coefficients')
plt.title('HenyeyGreenstein Fourier approximation')
plt.savefig("Graph.png")
