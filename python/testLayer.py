from layer import layer
from gaussLobatto import gausLobatto

import matplotlib.pyplot as plt
import layerlab as ll
import numpy as np

n = 64
mu, w = gausLobatto(n) #Valid


m = 12
eta = complex(1.1, 0.0)
alpha = 0.6
run = 0
if run == 0:
    diffuseLayerR = layer(mu, w, m)
    diffuseLayerR.setDiffuse(1.0)
    llLayerR = ll.Layer(mu, w, m)
    SMR = np.copy(diffuseLayerR.scatteringMatrix[:, :, 0])
    ll.setScatteringMatrix(llLayerR, SMR, 0)

    diffuseLayerG = layer(mu, w, m)
    diffuseLayerG.setDiffuse(0.0)
    llLayerG = ll.Layer(mu, w, m)
    SMG = np.copy(diffuseLayerG.scatteringMatrix[:, :, 0])
    ll.setScatteringMatrix(llLayerG, SMG, 0)

    diffuseLayerB = layer(mu, w, m)
    diffuseLayerB.setDiffuse(0.0)
    llLayerB = ll.Layer(mu, w, m)
    SMB = np.copy(diffuseLayerB.scatteringMatrix[:, :, 0])
    ll.setScatteringMatrix(llLayerB, SMB, 0)

    output = [llLayerR, llLayerG, llLayerB]
    print("Writing to disk..")
    storage = ll.BSDFStorage.fromLayerRGB("diffuse.bsdf", *output)
    storage.close()

elif run == 1:
    diffuseLayer = layer(mu, w, m)
    diffuseLayer.setDiffuse(0.8)
    temp = diffuseLayer.scatteringMatrix[:, :, 0]

    plt.figure()
    plt.imshow(temp)
    plt.savefig('images/diffuse' + str(0) + '.png')

    diffuseLayer = ll.Layer(mu, w, m)
    diffuseLayer.setDiffuse(0.8)
    i = 0
    topRow = np.concatenate((diffuseLayer[i].transmissionBottomTop, diffuseLayer[i].reflectionTop), axis=1)
    bottomRow = np.concatenate((diffuseLayer[i].reflectionBottom, diffuseLayer[i].transmissionTopBottom), axis=1)

    SM = np.concatenate((topRow, bottomRow), axis=0)
    print()

    plt.figure()
    plt.imshow(SM)
    plt.savefig('images/lldiffuse' + str(i) + '.png')
    plt.figure()
    plt.imshow(np.abs(SM - temp))
    plt.savefig('images/zerrdiffuse' + str(i) + '.png')
    print("error: ", str(np.sum(np.abs(SM - temp))))

elif run == 2:

    coating = layer(mu, w, m)
    print("setting microfacet")
    coating.setMicrofacet(eta, alpha)

    for i in range(m):
        plt.figure()
        temp = coating.scatteringMatrix[:, :, i]

        plt.imshow(temp)
        plt.savefig('images/microfacet' + str(i) + '.png')


    coatingll = ll.Layer(mu, w, m)
    coatingll.setMicrofacet(eta=eta, alpha=alpha, conserveEnergy=False)
    for i in range(m):
        topRow = np.concatenate((coatingll[i].transmissionBottomTop, coatingll[i].reflectionTop), axis=1)
        bottomRow = np.concatenate((coatingll[i].reflectionBottom, coatingll[i].transmissionTopBottom), axis=1)
        SM = np.concatenate((topRow, bottomRow), axis=0)


        plt.figure()
        plt.imshow(SM)
        plt.savefig('images/llmicrofacet' + str(i) + '.png')

    for i in range(m):
        topRow = np.concatenate((coatingll[i].transmissionBottomTop, coatingll[i].reflectionTop), axis=1)
        bottomRow = np.concatenate((coatingll[i].reflectionBottom, coatingll[i].transmissionTopBottom), axis=1)
        SM = np.concatenate((topRow, bottomRow), axis=0)
        temp = coating.scatteringMatrix[:, :, i]
        plt.figure()
        plt.imshow(np.abs(SM - temp))
        plt.savefig('images/errmicrofacet' + str(i) + '.png')

        # compute quartet errrors
        TtbErr = np.sum(np.abs(coating.getTtb(i) - coatingll[i].transmissionTopBottom))
        RbErr = np.sum(np.abs(coating.getRb(i) - coatingll[i].reflectionBottom))
        Rt = np.sum(np.abs(coating.getRt(i) - coatingll[i].reflectionTop))
        Tbt = np.sum(np.abs(coating.getTbt(i) - coatingll[i].transmissionBottomTop))
        globalErr = np.sum(np.abs(SM - temp))
        print("Global = " + str(globalErr))
        print (str(TtbErr) + ", " + str(RbErr))
        print (str(Rt) + ", " + str(TtbErr))
elif run == 3:

    coating = layer(mu, w, m)
    coating.setMicrofacet(eta, alpha)

    base = layer(mu, w, m)
    base.setDiffuse(0.8)

    layer.addToTop(base, coating)

elif run == 4:

    llLayerR = ll.Layer(mu, w, m)
    llLayerR.setDiffuse(1.0)
    
    llLayerG = ll.Layer(mu, w, m)
    llLayerG.setDiffuse(0.0)

    llLayerB = ll.Layer(mu, w, m)
    llLayerB.setDiffuse(0.0)
    output = [llLayerR, llLayerG, llLayerB]
    print("Writing to disk..")
    storage = ll.BSDFStorage.fromLayerRGB("diffusell.bsdf", *output)
    storage.close()
