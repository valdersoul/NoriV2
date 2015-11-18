/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/microfacet.h>

NORI_NAMESPACE_BEGIN

/// Ideal RoughConductor BRDF
class RoughDielectric : public BSDF {
public:
    RoughDielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);
    }


    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle) return Color3f(0.0f);
        /* Discrete BRDFs always evaluate to zero in Nori */
        float wiDotN = bRec.wi.z();
        float woDotN = bRec.wo.z();
        /* Determine the type of interaction */
        bool reflect = wiDotN*woDotN > 0.0f;
        float eta = Frame::cosTheta(bRec.wi) > 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;

        // the mirco scale normal
        Vector3f H;

        if (reflect) {
            /* Calculate the reflection half-vector */
            H = (bRec.wo+bRec.wi).normalized();
         } else {
            H = (bRec.wi + bRec.wo*eta).normalized();
        }

        /* Ensure that the half-vector points into the
           same hemisphere as the macrosurface normal */
        H *= sign(Frame::cosTheta(H));

        // Evaluate the microfacet normal distribution
        const float D = Microfacet::D(m_alpha, H);
        if (D== 0.0f) {
            return Color3f(0.0f);
        }

        float wiDotH = bRec.wi.dot(H);

        // Fresnel factor
        float cosThetaI;
        const float F = dielectricReflectance(wiDotH, cosThetaI, eta);

        /* Smith's shadow-masking function */
        const float G = Microfacet::G(m_alpha, bRec.wi, bRec.wo, H);


        if (reflect) {
            float fr = (F * G * D)/ (4.0f * std::abs(Frame::cosTheta(bRec.wi)));
            return Color3f(fr);
        } else {            


            /* Calculate the total amount of transmission */
            float sqrtDenom = bRec.wi.dot(H) + eta * bRec.wo.dot(H);
            float value = ((1 - F) * D * G * eta * eta
                * bRec.wi.dot(H) * bRec.wo.dot(H)) /
                (Frame::cosTheta(bRec.wi) * sqrtDenom * sqrtDenom);

            return std::abs(value);
        }
    }
    static inline float sign(float x) {
        if (x >= 0.0) return 1.0;
        else return -1.0;

    }

    static Vector3f refract(const Vector3f &wi, const Vector3f &n, float eta, float cosThetaT) {
        if (cosThetaT < 0)
            eta = 1 / eta;

        return n * (wi.dot(n) * eta + cosThetaT) - wi * eta;
    }



    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle)
            return 0.0f;
        /* Discrete BRDFs always evaluate to zero in Nori */

        bool reflect    = Frame::cosTheta(bRec.wi)
                        * Frame::cosTheta(bRec.wo) > 0.0f;

        float eta = Frame::cosTheta(bRec.wi) > 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        // the mirco scale normal
        Vector3f H;
        float dwh_dwo;

        if (reflect) {
            H    = (bRec.wi + bRec.wo).normalized();

            /* Jacobian of the half-direction mapping */
            dwh_dwo = 1.0f / (4.0f * bRec.wo.dot(H));
        } else {
            H    = (bRec.wi + bRec.wo * eta).normalized();

            /* Jacobian of the half-direction mapping */
            float sqrtDenom = bRec.wi.dot(H) + eta * bRec.wo.dot(H);
            dwh_dwo = (eta*eta * bRec.wo.dot(H)) / (sqrtDenom*sqrtDenom);
        }

        /* Ensure that the half-vector points into the
           same hemisphere as the macrosurface normal */
        H *= sign(Frame::cosTheta(H));

        float pm = Microfacet::pdf(m_alpha, H);

        float F = dielectricReflectance(bRec.wi.dot(H), m_intIOR / m_extIOR);
        pm *= reflect ? F : (1-F);

        return std::abs(pm * dwh_dwo);
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {        

        // sample the mirco scale normal
        Vector3f m = Microfacet::sample(m_alpha, sample);

        float pdf = Microfacet::pdf(m_alpha, m);
        if(pdf == 0.0) {
            std::cout << "sample pdf should nearly never be zero";
            return Color3f(0.0f);
        }

        float wiDotM = bRec.wi.dot(m);
        float cosThetaT = 0.0f;
        float F = dielectricReflectance(wiDotM, cosThetaT, m_intIOR / m_extIOR);


        bool sampleReflection = true;

        // check if reflection or refraction
        // for reuseable sample
        if(bRec.sampler->next1D() > F) {
            sampleReflection = false;
        }

        if (sampleReflection) {
            bRec.wo = 2.0f * wiDotM * m - bRec.wi;;
            bRec.eta = 1.0f;

            /* Side check */
            if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
                return Color3f(0.0f);

        } else {
            if (cosThetaT == 0)
                return Color3f(0.0f);
            bRec.wo = refract(bRec.wi, m, m_intIOR / m_extIOR, cosThetaT);
            bRec.eta = cosThetaT < 0 ? m_intIOR / m_extIOR : m_extIOR / m_extIOR;

            /* Side check */
            if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
                return Color3f(0.0f);


        }

        float G = Microfacet::G(m_alpha, bRec.wi, bRec.wo, m);
        float D = Microfacet::D(m_alpha, m);

        return std::abs(D * G * bRec.wi.dot(m) / (pdf * Frame::cosTheta(bRec.wi)));
    }

    std::string toString() const {
        return tfm::format(
            "roughDielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "  m_alpha = %f,\n"
            "]",
            m_intIOR, m_extIOR, m_alpha);
    }
    bool isDeltaBSDF() const {
        return false;
    }
private:
    float m_intIOR;
    float m_extIOR;

    // the roughness of the material
    float m_alpha;
};

NORI_REGISTER_CLASS(RoughDielectric, "roughDielectric");
NORI_NAMESPACE_END
