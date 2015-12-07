import numpy as np
import layerlab as ll
from microfacet import microfacetFourierSeries


class layer:
    def __init__(self, mu, w, m):
        self.nodes = mu
        self.weights = w
        self.modes = m
        self.scatteringMatrix = None

        nHalf = len(self.weights) / 2
        #
        self.weights[:nHalf] = (self.weights[:nHalf])[::-1]
        self.nodes[:nHalf] = (self.nodes[:nHalf])[::-1]

    def resolution(self):
        return len(self.nodes)

    def fourierOrders(self):
        return self.modes

    def getTtb(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[:resHalf, :resHalf, k]

    def getRb(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[:resHalf, resHalf:, k]

    def getRt(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[resHalf:, :resHalf, k]

    def getTbt(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[resHalf:, resHalf:, k]

    def setTtb(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[:resHalf, :resHalf, k] = M

    def setRb(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[:resHalf, resHalf:, k] = M

    def setRt(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[resHalf:, :resHalf, k] = M

    def setTbt(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[resHalf:, resHalf:, k] = M

    def scaleColumns(self, d, k):
        # first half
        scale = np.diag(d[:len(d) / 2])
        Tbt = self.getTbt(k).dot(scale)
        Rb = self.getRb(k).dot(scale)

        # second half
        scale = np.diag(d[len(d) / 2:])
        Ttb = self.getTtb(k).dot(scale)
        Rt = self.getRt(k).dot(scale)

        # set
        self.setRb(Rb, k)
        self.setTbt(Tbt, k)
        self.setTtb(Ttb, k)
        self.setRt(Rt, k)

    def applyWeights(self, F, l):
        WMarr = self.weights * np.abs(self.nodes)
        if l == 0:
            WMarr *= 2.0 * 3.14159265358979323846
        else:
            WMarr *= 3.14159265358979323846
        self.scatteringMatrix[:, :, l] = F
        self.scaleColumns(WMarr, l)
        return self.scatteringMatrix[:, :, l]
        # return F.dot(np.diag(temp))

    def setDiffuse(self, albedo):
        n = self.resolution()
        F = np.zeros((n, n))
        h = int(n / 2)
        self.scatteringMatrix = np.zeros((n, n, 1))
        for i in range(n):
            for o in range(n):
                if ((i < h and o >= h) or (o < h and i >= h)):
                    F[o, i] = albedo * 0.31830988618379067154

        self.scatteringMatrix[:, :, 0] = self.applyWeights(F, 0)

    def setMicrofacet(self, eta, alpha):
        n = self.resolution()
        h = int(n / 2)
        fourierOrdersTarget = self.fourierOrders()
        FL = np.zeros((n, n, fourierOrdersTarget))
        for i in range(n):
            for o in range(n):
                coeffs = ll.microfacetFourierSeries(self.nodes[o], self.nodes[i], eta, alpha, fourierOrdersTarget,
                                                    10e-4)
                # coeffs = microfacetFourierSeries( self.nodes[o],  self.nodes[i], eta, alpha, fourierOrdersTarget)
                for l in range(len(coeffs[:fourierOrdersTarget])):
                    FL[o, i, l] = coeffs[l]

        self.scatteringMatrix = np.zeros((n, n, fourierOrdersTarget))
        for i in range(fourierOrdersTarget):
            self.scatteringMatrix[:, :, i] = self.applyWeights(FL[:, :, i], i)
