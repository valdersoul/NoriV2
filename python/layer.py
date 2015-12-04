import numpy as np
import layerlab as ll
from microfacet import microfacetFourierSeries


class layer:
    def __init__(self, mu, w, m):
        self.nodes = mu
        self.weights = w
        self.modes = m
        self.scatteringMatrix = None

    def resolution(self):
        return len(self.nodes)

    def fourierOrders(self):
        return self.modes

    def applyWeights(self, F, l):
        temp = self.weights * np.abs(self.nodes)
        if l == 0:
            temp *= 2.0 * np.pi
        else:
            temp *= np.pi
        return F.dot(np.diag(temp))

    def setDiffuse(self, albedo):
        n = self.resolution()
        F = np.zeros((n, n))
        h = int(n / 2)
        self.scatteringMatrix = np.zeros((n, n, 1))
        for i in range(n):
            for o in range(n):
                if (i < h <= o) or (o < h <= i):
                    F[o, i] = albedo / np.pi

        self.scatteringMatrix[:, :, 0] = self.applyWeights(F, 0)

    def setMicrofacet(self, eta, alpha):
        n = self.resolution()
        h = int(n / 2)
        fourierOrdersTarget = self.fourierOrders()
        FL = np.zeros((n, n, fourierOrdersTarget))
        for i in range(n):
            for o in range(n):
                coeffs = ll.microfacetFourierSeries(- self.nodes[o], - self.nodes[i], eta, alpha, fourierOrdersTarget, 10e-4)
                #coeffs = microfacetFourierSeries(- self.nodes[o], - self.nodes[i], eta, alpha, fourierOrdersTarget)
                for l in range(len(coeffs[:fourierOrdersTarget])):
                    FL[o, i, l] = coeffs[l]

        SM = np.zeros((n, n, fourierOrdersTarget))
        for i in range(fourierOrdersTarget):
            SM[:, :, i] = self.applyWeights(FL[:, :, i], i)

        self.scatteringMatrix = SM
