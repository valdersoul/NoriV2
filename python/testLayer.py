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
run = 2
if run == 1:
    diffuseLayer = layer(mu, w, m)
    diffuseLayer.setDiffuse(0.8)
    temp = diffuseLayer.scatteringMatrix[:, :, 0]
    resHalf = n / 2
    #temp[resHalf:, :resHalf] = np.fliplr(temp[resHalf:, :resHalf])
    #temp[:resHalf, resHalf:] = np.flipud(temp[:resHalf, resHalf:])
    plt.figure()
    plt.imshow(temp)
    plt.savefig('images/diffuse' + str(0) + '.png')

    diffuseLayer = ll.Layer(mu, w, m)
    diffuseLayer.setDiffuse(0.8)
    i = 0
    topRow = np.concatenate((diffuseLayer[i].transmissionTopBottom, diffuseLayer[i].reflectionBottom), axis=1)
    bottomRow = np.concatenate((diffuseLayer[i].reflectionTop, diffuseLayer[i].transmissionBottomTop), axis=1)
    SM = np.concatenate((topRow, bottomRow), axis=0)

    plt.figure()
    plt.imshow(SM)
    plt.savefig('images/lldiffuse' + str(i) + '.png')
    plt.figure()
    plt.imshow(np.abs(SM - temp))
    plt.savefig('images/zerrdiffuse' + str(i) + '.png')
    print(np.sum(np.abs(SM - temp)))

elif run == 2:

    coating = layer(mu, w, m)
    print("setting microfacet")
    coating.setMicrofacet(eta, alpha)

    for i in range(m):
        print(i)
        plt.figure()
        temp = coating.scatteringMatrix[:, :, i]
        resHalf = n / 2
        #temp[resHalf:, :resHalf] = np.fliplr(temp[resHalf:, :resHalf])
        #temp[:resHalf, resHalf:] = np.flipud(temp[:resHalf, resHalf:])
        plt.imshow(temp)
        plt.savefig('images/microfacet' + str(i) + '.png')


    coatingll = ll.Layer(mu, w, m)
    coatingll.setMicrofacet(eta=eta, alpha=alpha, conserveEnergy=False)
    for i in range(m):
        SM = coatingll.matrix(i)


        plt.figure()
        plt.imshow(SM)
        plt.savefig('images/llmicrofacet' + str(i) + '.png')

    for i in range(m):
        topRow = np.concatenate((coatingll[i].transmissionTopBottom, coatingll[i].reflectionBottom), axis=1)
        bottomRow = np.concatenate((coatingll[i].reflectionTop, coatingll[i].transmissionBottomTop), axis=1)
        SM = np.concatenate((topRow, bottomRow), axis=0)
        #SM = coatingll.matrix(i)
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

