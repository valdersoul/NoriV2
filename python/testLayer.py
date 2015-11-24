from layer import layer
from gaussLobatto import gausLobatto

import matplotlib.pyplot as plt
import layerlab as ll
import numpy as np

n = 128
mu, w = gausLobatto(n)

m = 20
eta = complex(1.1, 0.0)
alpha = 0.6
run = 2
if run == 1:
    diffuseLayer = layer(mu, w, m)
    diffuseLayer.setDiffuse(0.8)
    plt.figure()
    plt.imshow(diffuseLayer.scatteringMatrix[:, :, 0])
    plt.savefig('images/diffuse' + str(0) + '.png')

    diffuseLayer = ll.Layer(mu, w, m)
    diffuseLayer.setDiffuse(0.8)
    i = 0
    topRow = np.concatenate((diffuseLayer[i].transmissionBottomTop, diffuseLayer[i].reflectionTop), axis=1)
    bottomRow = np.concatenate((diffuseLayer[i].reflectionBottom, diffuseLayer[i].transmissionTopBottom), axis=1)
    SM = np.concatenate((topRow, bottomRow), axis=0)

    plt.figure()
    plt.imshow(SM)
    plt.savefig('images/lldiffuse' + str(i) + '.png')

elif run == 2:

    coating = layer(mu, w, m)
    coating.setMicrofacet(eta, alpha)

    for i in range(m):
        plt.figure()
        plt.imshow(coating.scatteringMatrix[:, :, i], cmap=plt.get_cmap('gray'))
        plt.savefig('images/microfacet' + str(i) + '.png')

    coating = ll.Layer(mu, w, m)
    coating.setMicrofacet(eta=eta, alpha=alpha)
    for i in range(m):
        topRow = np.concatenate((coating[i].transmissionBottomTop, coating[i].reflectionTop), axis=1)
        bottomRow = np.concatenate((coating[i].reflectionBottom, coating[i].transmissionTopBottom), axis=1)
        SM = np.concatenate((topRow, bottomRow), axis=0)

        plt.figure()
        plt.imshow(SM, cmap=plt.get_cmap('gray'))
        plt.savefig('images/llmicrofacet' + str(i) + '.png')
