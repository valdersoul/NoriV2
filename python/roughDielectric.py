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
eta_bot = 1.5
alpha_bot = 0.2  # Beckmann roughness of bottom layer (gold)

# Construct quadrature scheme suitable for the material
n, m = ll.parameterHeuristicMicrofacet(eta=eta_bot, alpha=alpha_bot)
#n = 64
#m = 200
mu, w = ll.quad.gaussLobatto(n)
print("# of nodes = %i, fourier orders = %i" % (n, m))

lloutput = []

# Construct diffuse bottom layer for each channel
print("2. Creating metal layer")
lR = ll.Layer(mu, w, m)
lR.setMicrofacet(eta=eta_bot, alpha=alpha_bot)
lG = ll.Layer(mu, w, m)
lG.setMicrofacet(eta=eta_bot, alpha=alpha_bot)
lB = ll.Layer(mu, w, m)
lB.setMicrofacet(eta=eta_bot, alpha=alpha_bot)
lloutput= [lR, lG, lB]

# .. and write to disk
print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("rDielectricLL.bsdf", *lloutput )
storage.close()


#####################################################################

output = []

# Construct diffuse bottom layer for each channel
print("1. Creating metal layer")
l = layer(mu, w, m)
l.setMicrofacet(eta_bot, alpha_bot)

# Apply coating
print("1. Applying coating..")
output = [l, l, l]

outputForReal = []
for channel in range(3):
    llLayer = ll.Layer(mu, w, m)
    for i in range(m):
        SM = np.copy(output[channel].scatteringMatrix[:, :, i])

        ll.setScatteringMatrix(llLayer, SM, i)
    outputForReal.append(llLayer)
    

print("Writing to disk..")
storage = ll.BSDFStorage.fromLayerRGB("rDielectric.bsdf", *outputForReal)
storage.close()





























