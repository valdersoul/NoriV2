from layer import layer
from gaussLobatto import gausLobatto

import matplotlib.pyplot as plt
import layerlab as ll
import numpy as np
from utils.materials import gold
from utils.cie import get_rgb

import sys
sys.path.append('.')

def getMatrix(output, i):
    topRow = np.concatenate((output[i].transmissionBottomTop, output[i].reflectionTop), axis=1)
    bottomRow = np.concatenate((output[i].reflectionBottom, output[i].transmissionTopBottom), axis=1)
    SM = np.concatenate((topRow, bottomRow), axis=0)
    return SM

def printMatrix(SM1, SM2, i):
    plt.figure()
    plt.subplot(1,2,1)
    plt.title('OUR')
    plt.imshow(SM1)
    plt.subplot(1,2,2)
    plt.title('ll')
    plt.imshow(SM2)
    plt.savefig('images/addComp' + str(i) + '.png')

albedo = 0.5
eta_top = 1.5
alpha  = 0.02  # Beckmann roughness

# Construct quadrature scheme suitable for the material
print('Computing gold IOR parameters')
eta_bot = get_rgb(gold)

alpha_top = 0.1  # Beckmann roughness of top layer (coating)
alpha_bot = 0.1  # Beckmann roughness of bottom layer (gold)

# Construct quadrature scheme suitable for the material
n_top, m_top = ll.parameterHeuristicMicrofacet(eta=eta_top, alpha=alpha_top)
n_bot, m_bot = ll.parameterHeuristicMicrofacet(eta=eta_bot[0], alpha=alpha_bot)
n = max(n_top, n_bot)  # Max of zenith angle discretization
m = m_top              # Number of Fourier orders determined by top layer
mu, w = ll.quad.gaussLobatto(n)
print("# of nodes = %i, fourier orders = %i" % (n, m))


print("1. Creating coating layer")
coating = layer(mu, w, m)
coating.setMicrofacet(eta_top, alpha_top)


output = []
for channel in range(3):
    # Construct diffuse bottom layer for each channel
    print("1. Creating metal layer")
    l = layer(mu, w, m)
    l.setMicrofacet(eta_bot[channel], alpha_bot)

    # Apply coating
    print("1. Applying coating..")
    output.append(layer.addToTop(coating, l))

outputForReal = []
for channel in range(3):
    llLayer = ll.Layer(mu, w, m)
    for i in range(m):
        SM = np.copy(output[channel].scatteringMatrix[:, :, i])

        ll.setScatteringMatrix(llLayer, SM, i)
    outputForReal.append(llLayer)

storage = ll.BSDFStorage.fromLayerRGB("ourOutput.bsdf", *outputForReal)
storage.close()


# Construct coating layer
print("2. Creating coating layer")
coating = ll.Layer(mu, w, m)
coating.setMicrofacet(eta=eta_top, alpha=alpha_top)

lloutput = []
for channel in range(3):
    # Construct diffuse bottom layer for each channel
    print("2. Creating metal layer")
    l = ll.Layer(mu, w, m)
    l.setMicrofacet(eta=eta_bot[channel], alpha=alpha_bot)

    # Apply coating
    print("2. Applying coating..")
    l.addToTop(coating)
    lloutput .append(l)

# .. and write to disk
print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("lloutput.bsdf", *lloutput )
storage.close()


for i in range(12):
    SM1 = output[0].scatteringMatrix[:, :, i]
    SM2 = getMatrix(lloutput[0], i)
    printMatrix(SM1, SM2, i)


