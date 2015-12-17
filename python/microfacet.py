import frame as f
import math
import numpy as np
import layerlab as ll
from BSDFFourier import expcos_fseries, fseries_convolve, expcosCoefficientCount
from filon import filonIntegrate
import scipy.special as sp
import sys
from scipy.linalg import svd



def fresnelDielectric(cosThetaI_, cosThetaT_, eta):
    if eta == 1:
        cosThetaT_ = -cosThetaI_
        return 0.0, cosThetaT_

    # /* Using Snell's law, calculate the squared sine of the angle between the normal and the transmitted ray */
    scale = eta
    if cosThetaI_ > 0:
        scale = 1 / eta
    cosThetaTSqr = 1 - (1 - cosThetaI_ * cosThetaI_) * (scale * scale)

    # /* Check for total internal reflection */
    if cosThetaTSqr <= 0.0:
        cosThetaT_ = 0.0
        return 1.0, cosThetaT_

    # /* Find the absolute cosines of the incident/transmitted rays */
    cosThetaI = np.abs(cosThetaI_)
    cosThetaT = math.sqrt(cosThetaTSqr)

    Rs = (cosThetaI - eta * cosThetaT) / (cosThetaI + eta * cosThetaT)
    Rp = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT)

    cosThetaT_ = cosThetaT
    if cosThetaI_ > 0:
        cosThetaT_ = -cosThetaT

    # /* No polarization -- return the unpolarized reflectance */
    return 0.5 * (Rs * Rs + Rp * Rp), cosThetaT_


def fresnelConductor(cosThetaI, eta_):
    eta = eta_.real
    k = eta_.imag

    cosThetaI2 = cosThetaI * cosThetaI
    sinThetaI2 = 1.0 - cosThetaI2
    sinThetaI4 = sinThetaI2 * sinThetaI2

    temp1 = eta * eta - k * k - sinThetaI2
    a2pb2 = f.safe_sqrt(temp1 * temp1 + 4 * k * k * eta * eta)
    a = f.safe_sqrt(0.5 * (a2pb2 + temp1))

    term1 = a2pb2 + cosThetaI2
    term2 = 2 * a * cosThetaI

    Rs2 = (term1 - term2) / (term1 + term2)

    term3 = a2pb2 * cosThetaI2 + sinThetaI4
    term4 = term2 * sinThetaI2

    Rp2 = Rs2 * (term3 - term4) / (term3 + term4)

    return 0.5 * (Rp2 + Rs2)


def smithG1(v, m, alpha):
    tanTheta = np.abs(f.tanTheta(v))

    # /* Can't see the back side from the front and vice versa */
    if np.dot(v, m) * f.cosTheta(v) <= 0:
        return 0.0
    # if alpha * tanTheta == 0.0:
    #    return 0.0
    a = 1.0 / (alpha * tanTheta)
    if a < 1.6:
        aSqr = a * a
        return (3.535 * a + 2.181 * aSqr) / (1.0 + 2.276 * a + 2.577 * aSqr)

    return 1.0


def getAB(mu_i, mu_o, eta, alpha):
    reflect = (-mu_i * mu_o) > 0
    sinMu2 = f.safe_sqrt((1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o))

    if reflect:
        temp = 1.0 / (alpha * (mu_i - mu_o))
        A = (mu_i * mu_i + mu_o * mu_o - 2.0) * temp * temp
        B = 2.0 * sinMu2 * temp * temp
    else:
        temp = 1.0 / (alpha * (mu_i - eta.real * mu_o))
        A = (mu_i * mu_i - 1.0 + eta.real * eta.real * (mu_o * mu_o - 1.0)) * temp * temp
        B = 2 * eta.real * sinMu2 * temp * temp

    return A, B


# Float mu_o
# Float mu_i
# std::complex<Float> eta_ .. here cmath.complex
# Float alpha
# Float phi_d
def microfacetNoExp(mu_o, mu_i, eta_, alpha, phi_d):
    sinThetaI = f.safe_sqrt(1.0 - mu_i * mu_i)
    sinThetaO = f.safe_sqrt(1.0 - mu_o * mu_o)
    cosPhi = math.cos(phi_d)
    sinPhi = math.sin(phi_d)

    wi = np.array([-sinThetaI, 0.0, -mu_i])
    wo = np.array([sinThetaO * cosPhi, sinThetaO * sinPhi, mu_o])
    reflect = (-mu_i * mu_o) > 0

    if (mu_o == 0.0) or (mu_i == 0.0):
        return 0.0

    conductor = eta_.imag != 0.0
    if conductor and not reflect:
        return 0.0

    eta = eta_
    if not (-mu_i > 0 or conductor):
        eta = (complex(1.0, 0.0) / eta_)

    scalingFactor = 1.0
    if not reflect:
        scalingFactor = eta.real
    H = wi + wo * scalingFactor
    H /= np.linalg.norm(H)

    H *= np.sign(H[2])  # math::signum(Frame::cosTheta(H));

    cosThetaH2 = H[2] * H[2]
    D = 1.0 / (math.pi * alpha * alpha * cosThetaH2 * cosThetaH2)

    # Calculate the fresnel term
    F = 0.0
    if not conductor:
        F, _ = fresnelDielectric(np.dot(wi, H), 0, eta_.real)
    else:
        F = fresnelConductor(np.abs(np.dot(wi, H)), eta)
    G = smithG1(wi, H, alpha) * smithG1(wo, H, alpha)

    if reflect:
        return F * D * G / (4.0 * np.abs(mu_i * mu_o))
    else:
        sqrtDenom = np.dot(wi, H) + eta.real * np.dot(wo, H)

        return np.abs(((1 - F) * D * G * eta.real * eta.real * np.dot(wi, H) * np.dot(wo, H)) / (
            mu_i * mu_o * sqrtDenom * sqrtDenom))


def Bmax(n, relerr):
    if relerr >= 1e-1:
        return 0.1662 * math.pow(n, 2.05039)
    elif relerr >= 1e-2:
        return 0.0818 * math.pow(n, 2.04982)
    elif relerr >= 1e-3:
        return 0.0538 * math.pow(n, 2.05001)
    elif relerr >= 1e-4:
        return 0.0406 * math.pow(n, 2.04686)
    elif relerr >= 1e-5:
        return 0.0337 * math.pow(n, 2.03865)
    elif relerr >= 1e-6:
        return 0.0299 * math.pow(n, 2.02628)
    else:
        print("Bmax(): unknown relative error bound!")


def microfacetNoExpFourierSeries(mu_o, mu_i, etaC, alpha, n, phiMax):
    # const
    nEvals = 200
    eps = sys.float_info.epsilon

    def foo(a):
        return microfacetNoExp(mu_o, mu_i, etaC, alpha, a)

    reflect = - mu_i * mu_o > 0.0
    sinMu2 = f.safe_sqrt((1.0 - mu_i * mu_i) * (1.0 - mu_o * mu_o))
    phiCritical = 0.0
    conductor = etaC.imag != 0.0

    if -mu_i > 0.0 or conductor:
        eta = etaC
    else:
        eta = complex(1.0, 0.0) / etaC

    if reflect:
        if not conductor:
            if sinMu2 == 0:
                print("reflect sinMu2 == 0")
                temp = -1.0
            else:
                temp = (2.0 * eta.real * eta.real - mu_i * mu_o - 1.0) / sinMu2
            phiCritical = f.safe_acos(temp)
    elif not reflect:
        if eta.real > 1.0:
            etaDenser = eta.real
        else:
            etaDenser = 1.0 / eta.real
        if sinMu2 == 0:
            print("refract sinMu2 == 0")
            temp = -1.0
        else:
            temp = (1.0 - etaDenser * mu_i * mu_o) / (etaDenser * sinMu2)
        phiCritical = f.safe_acos(temp)
    if not conductor and phiCritical > eps and phiCritical < np.pi - eps and phiCritical < phiMax - eps:
        #   Uh oh, some high frequency content leaked in the generally low frequency part.
        #   Increase the number of coefficients so that we can capture it. Fortunately, this
        #   happens very rarely.
        print("Uh oh case")
        n = max(n, 100)
    # validated and tested length and value
    if reflect:
        if phiCritical > eps and phiCritical < phiMax - eps:

            b1 = filonIntegrate(foo, n, nEvals, 0.0, phiCritical)
            b2 = filonIntegrate(foo, n, nEvals, phiCritical, phiMax)
            coeffs = np.concatenate((b1, b2), axis=0)
        else:
            coeffs = filonIntegrate(foo, n, nEvals, 0.0, phiMax)
    else:
        coeffs = filonIntegrate(foo, n, nEvals, 0.0, min(phiCritical, phiMax))

    if phiMax < np.pi - eps:
        #  /* Precompute some sines and cosines */
        cosPhi = np.zeros(n)
        sinPhi = np.zeros(n)
        for i in range(n):
            sinPhi[i] = np.sin(i * phiMax)
            cosPhi[i] = np.cos(i * phiMax)
        # tested and validated values over sum
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1):
                if i != j:
                    A[i, j] = A[j, i] = (i * cosPhi[j] * sinPhi[i] - j * cosPhi[i] * sinPhi[j]) / (i * i - j * j)
                elif i != 0:
                    A[i, i] = (np.sin(2.0 * i * phiMax) + 2.0 * i * phiMax) / (4.0 * i)
                else:
                    A[i, i] = phiMax

        U, sigma, V = svd(A)
        V = V.T

        if sigma[0] == 0:
            return [0.0]

        temp = np.zeros(n)
        coeffs[0] *= np.pi
        coeffs[1:n] *= 0.5 * np.pi
        for i in range(n):
            if sigma[i] < 1e-9 * sigma[0]:
                break
            if len(coeffs) == len(U[:, i]):
                temp += V[:, i] * np.dot(U[:, i], coeffs) / sigma[i]
            else:
                print("Error shapes doesnt match")
                print("len(coeffs)" + str(len(coeffs)))
                print("len(U[:, i])" + str(len(U[:, i])))
        coeffs = temp

    return coeffs


def microfacetFourierSeries(mu_o, mu_i, etaC, alpha, n, relerr=10e-4):
    conductor = (etaC.imag != 0.0)
    reflect = - mu_o * mu_i > 0.0

    # numerical errors
    if conductor and not reflect:
        return [0.0]

    # validated visual
    A, B = getAB(mu_i, mu_o, etaC, alpha)

    # validate and tested
    if sp.i0e(B) * np.exp(A + B) < 1e-10:
        return [0.0]

    # validated visual
    B_max = Bmax(n, relerr)
    if B > B_max:
        A = A + B - B_max + math.log(sp.i0e(B) / sp.i0e(B_max))
        B = B_max

    # validated and tested with length and value
    expcos_coeffs = expcos_fseries(A, B, relerr)

    if B == 0.0:
        phiMax = f.safe_acos(-1.0)
    else:
        phiMax = f.safe_acos(1.0 + np.log(relerr) / B)

    # validated length value is 0.00264447443559 PROBELM IS THE SVD computation

    lowfreq_coeffs = microfacetNoExpFourierSeries(mu_o, mu_i, etaC, alpha, 12, phiMax)

    # valid?
    result = fseries_convolve(lowfreq_coeffs, len(lowfreq_coeffs), expcos_coeffs, len(expcos_coeffs))

    for i in range(len(result)):
        if result[i] == 0 or np.abs(result[i]) < result[0] * relerr:
            result = result[:i]
            break

    return result
