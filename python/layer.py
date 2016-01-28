import numpy as np
import layerlab as ll
from microfacet import microfacetFourierSeries
import math

class layer:
    def __init__(self, mu, w, m):
        self.nodes = mu
        self.weights = w
        self.modes = m
        self.scatteringMatrix = None

        if self.nodes[0] < self.nodes[1]:
            nHalf = len(self.weights) / 2
            #
            self.weights[:nHalf] = (self.weights[:nHalf])[::-1]
            self.nodes[:nHalf] = (self.nodes[:nHalf])[::-1]

    @staticmethod
    def addToTop(l1, l2):
        n = l1.resolution() / 2
        I = np.identity(n)

        for i in range(l1.fourierOrders()):
            G_tb = (I - np.dot(l1.getRb(i), l2.getRt(i)))

            G_bt = (I - np.dot(l2.getRt(i), l1.getRb(i)))

            result = np.linalg.solve(G_tb, l1.getTtb(i))
            Ttb = np.dot(l2.getTtb(i), result)

            temp = np.dot(l1.getRb(i), l2.getTbt(i))
            result = np.linalg.solve(G_tb, temp)
            Rb = l2.getRb(i) + np.dot(l2.getTtb(i), result)

            result = np.linalg.solve(G_bt, l2.getTbt(i))
            Tbt = np.dot(l1.getTbt(i), result)

            temp = np.dot(l2.getRt(i), l1.getTtb(i))
            result = np.linalg.solve(G_bt, temp)
            Rt = l1.getRt(i) + np.dot(l1.getTbt(i), result)

            l2.setRt(Rt, i)
            l2.setRb(Rb, i)
            l2.setTtb(Ttb, i)
            l2.setTbt(Tbt, i)

        return l2


    def resolution(self):
        return len(self.nodes)

    def setScatteringMatrix(self, M):
        n = self.resolution()
        self.scatteringMatrix = np.zeros((n, n, self.fourierOrders()))
        for l in range(self.fourierOrders()):
            self.setTbt(M[l].transmissionBottomTop, l)
            self.setTtb(M[l].transmissionTopBottom, l)
            self.setRb(M[l].reflectionBottom, l)
            self.setRt(M[l].reflectionTop, l)

    def fourierOrders(self):
        return self.modes

    def getTbt(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[:resHalf, :resHalf, k]

    def getRt(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[:resHalf, resHalf:, k]

    def getRb(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[resHalf:, :resHalf, k]

    def getTtb(self, k):
        resHalf = self.resolution() / 2
        return self.scatteringMatrix[resHalf:, resHalf:, k]

    def setTbt(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[:resHalf, :resHalf, k] = M

    def setRt(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[:resHalf, resHalf:, k] = M

    def setRb(self, M, k):
        resHalf = self.resolution() / 2
        self.scatteringMatrix[resHalf:, :resHalf, k] = M

    def setTtb(self, M, k):
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
        self.scatteringMatrix = np.zeros((n, n, self.fourierOrders()))
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
                #coeffs = ll.microfacetFourierSeries(-self.nodes[o], -self.nodes[i], eta, alpha, fourierOrdersTarget, 1e-3) # works good but paper is different
                coeffs = microfacetFourierSeries(-self.nodes[o], -self.nodes[i], eta, alpha, fourierOrdersTarget, 1e-3)
                """
                if len(llcoeffs) == len(coeffs):
                    err = np.sum(np.abs(np.array(coeffs) - np.array(llcoeffs)))
                    if err > 1e-6:
                        print("absolute error = " + str(err) + " mu_o = " + str(self.nodes[o]) + " mu_i = " + str(self.nodes[i]))
                        print("======================================")
                else:
                    print("len(expcos_coeffs) = " + str(len(coeffs)))
                    print("len(llexpcos_coeffs) = " + str(len(llcoeffs)))
                    print(" mu_o = " + str(self.nodes[o]) + " mu_i = " + str(self.nodes[i]))
                    print("======================================")
                """

                for l in range(min(fourierOrdersTarget, len(coeffs))):
                    if math.isnan(coeffs[l]):
                        print("NAN found l = " + str(l) + " mu_o = " + str(self.nodes[o]) + " mu_i = " + str(self.nodes[i]))
                    FL[o, i, l] = coeffs[l]

        self.scatteringMatrix = np.zeros((n, n, fourierOrdersTarget))
        for i in range(fourierOrdersTarget):
            self.scatteringMatrix[:, :, i] = self.applyWeights(FL[:, :, i], i)
