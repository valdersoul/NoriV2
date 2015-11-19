import frame as f
import math
import numpy as np
import layerlab as ll


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
        A = (mu_i * mu_i + mu_o * mu_o - 2) * temp * temp
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
