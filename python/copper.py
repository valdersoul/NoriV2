from layer import layer
from gaussLobatto import gausLobatto

import matplotlib.pyplot as plt
import layerlab as ll
import numpy as np
from utils.materials import copper
from utils.cie import get_rgb

import sys
sys.path.append('.')


# Construct quadrature scheme suitable for the material
print('Computing copper IOR parameters')
eta_bot = get_rgb(copper)
alpha_bot = 0.1  # Beckmann roughness of bottom layer (gold)

# Construct quadrature scheme suitable for the material
n, m = ll.parameterHeuristicMicrofacet(eta=eta_bot[0], alpha=alpha_bot)
mu, w = ll.quad.gaussLobatto(n)
print("# of nodes = %i, fourier orders = %i" % (n, m))

lloutput = []
for channel in range(3):
    # Construct diffuse bottom layer for each channel
    print("2. Creating metal layer")
    l = ll.Layer(mu, w, m)
    l.setMicrofacet(eta=eta_bot[channel], alpha=alpha_bot)
    lloutput.append(l)

# .. and write to disk
print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("copperLL.bsdf", *lloutput )
storage.close()


#####################################################################

output = []
for channel in range(3):
    # Construct diffuse bottom layer for each channel
    print("1. Creating metal layer")
    l = layer(mu, w, m)
    l.setMicrofacet(eta_bot[channel], alpha_bot)

    # Apply coating
    print("1. Applying coating..")
    output.append(l)

outputForReal = []
for channel in range(3):
    llLayer = ll.Layer(mu, w, m)
    for i in range(m):
        SM = np.copy(output[channel].scatteringMatrix[:, :, i])

        ll.setScatteringMatrix(llLayer, SM, i)
    outputForReal.append(llLayer)
    

print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("ourOutputCopper.bsdf", *outputForReal)
storage.close()





























