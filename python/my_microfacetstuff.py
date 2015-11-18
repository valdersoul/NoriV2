import numpy as np

def safe_sqrt(value):
    return math.sqrt(math.max(0.0, value))

# Float mu_o
# Float mu_i
# std::complex<Float> eta_ .. here cmath.complex
# Float alpha
# Float phi_d
def microfacetNoExp(mu_o, mu_i, eta_, alpha, phi_d):
    sinThetaI = safe_sqrt(1-mu_i*mu_i)
    sinThetaO = safe_sqrt(1-mu_o*mu_o)
    cosPhi = math.cos(phi_d)
    sinPhi = math.sin(phi_d)

    wi = np.array([-sinThetaI, 0.0, -mu_i])
    wo = np.array([sinThetaO*cosPhi, sinThetaO*sinPhi, mu_o])
    reflect = (-mu_i*mu_o) > 0;

    if ((mu_o == 0) or (mu_i == 0)):
        return 0.0

    conductor = eta_.imag != 0.0
    if (conductor and not reflect):
        return 0.0

    eta = eta_
    if not (-mu_i > 0 or conductor):
        eta = (complex(1.0, 0.0) / eta_)        

    scalingFactor = 1.0
    if reflect:
        scalingFactor = eta.real
    H = wi + wo * scalingFactor
    H = H / norm(H)

    H = H * np.sign(H[2]) # math::signum(Frame::cosTheta(H));

    cosThetaH2 = H[2] * H[2]
    D = 1.0 / (math.pi * alpha * alpha * cosThetaH2 * cosThetaH2)

    # Calculate the fresnel term
    F = 0.0
    if not conductor:
        F = fresnelDielectric(np.dot(wi, H), eta_.real)
    else:
        F = fresnelConductor(math.abs(dot(wi, H)), eta)
    G = smithG1(wi, H, alpha) * smithG1(wo, H, alpha)

    if reflect:
        return F * D * G / (4.0 * math.abs(mu_i*mu_o))
    else:
        sqrtDenom = dot(wi, H) + eta.real * dot(wo, H)

        return math.abs(((1 - F) * D * G * eta.real * eta.real * dot(wi, H)
            * dot(wo, H)) / (mu_i*mu_o * sqrtDenom * sqrtDenom));