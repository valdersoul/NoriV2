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
        /* Discrete BRDFs always evaluate to zero in Nori */
        float wiDotN = bRec.wi.z();
        float woDotN = bRec.wo.z();

        bool reflect = wiDotN*woDotN >= 0.0f;

        float eta = wiDotN < 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        // the mirco scale normal
        Vector3f m;

        if (reflect)
            m = sign(wiDotN)*(bRec.wi + bRec.wo).normalized();
        else
            m = -(bRec.wi*eta + bRec.wo).normalized();

        float wiDotM = bRec.wi.dot(m);
        float woDotM = bRec.wo.dot(m);
        float cosThetaI;
        float F = dielectricReflectance(eta, wiDotM, cosThetaI);
        float G = Microfacet::G(m_alpha, bRec.wi, bRec.wo, m);
        float D = Microfacet::D(m_alpha, m);

        if (reflect) {
            float fr = (F*G*D*0.25f)/std::abs(wiDotN);
            return Color3f(fr);
        } else {
            float nomSquared = (eta*wiDotM + woDotM) * (eta*wiDotM + woDotM);
            float fs = std::abs(wiDotM*woDotM)*(1.0f - F)*G*D/(nomSquared*std::abs(wiDotN));
            return Color3f(fs);
        }

    }
    static float sign(float x) {
        if (x > 0.0) return 1.0;
        if (x < 0.0) return -1.0;
        return 0.0;
    }



    float pdf(const BSDFQueryRecord &bRec) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        float wiDotN = bRec.wi.z();
        float woDotN = bRec.wo.z();

        bool reflect = wiDotN*woDotN >= 0.0f;

        float eta = wiDotN < 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        // the mirco scale normal
        Vector3f m;

        if (reflect) {
            m = sign(wiDotN)*(bRec.wi + bRec.wo).normalized();
        } else {
            m = -(bRec.wi*eta + bRec.wo).normalized();
        }

        float wiDotM = bRec.wi.dot(m);
        float woDotM = bRec.wo.dot(m);

        float cosThetaT;
        float F = dielectricReflectance(eta, wiDotM, cosThetaT);
        float pm = Microfacet::pdf(m_alpha, m);


        float pdf;
        if (reflect) {
            pdf = F * pm * 0.25f / std::abs(wiDotM);
         } else {
            float denomSquared = (eta*wiDotM + woDotM) * (eta*wiDotM + woDotM);
            pdf = (1.0f - F) * pm * std::abs(woDotM)/denomSquared ;
        }
        return pdf;

    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {        
        float wiDotN = bRec.wi.z();

        float eta = wiDotN < 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;

        // sample the mirco scale normal
        Vector3f m = Microfacet::sample(m_alpha, sample);

        float pdf = Microfacet::pdf(m_alpha, m);
        if(pdf == 0.0) return Color3f(0.0f);

        float wiDotM = bRec.wi.dot(m);

        float cosThetaT = 0.0f;

        float F = dielectricReflectance(eta, wiDotM, cosThetaT);
        float etaM = wiDotM < 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;


        bool reflect = false;

        // check if reflection or refraction
        // for reuseable sample
        if(sample(0) < F) {
            reflect = true;
        }

        if (reflect) {
            bRec.wo = 2.0f * wiDotM * m - bRec.wi;
        } else {
            bRec.wo = (etaM * wiDotM - sign(wiDotM) * cosThetaT) * m - etaM * bRec.wi;
        }

        float woDotN = bRec.wo.z();

        bool reflected = wiDotN * woDotN > 0.0f;
        if (reflected != reflect) return Color3f(0.0f);

        float G = Microfacet::G(m_alpha, bRec.wi, bRec.wo, m);
        float D = Microfacet::D(m_alpha, m);


        Color3f col = Color3f(std::abs(wiDotM) * G * D/(std::abs(wiDotN) * pdf));
        if (reflect) {
            return col * F;
        }
        else
            return col * (1.0f - F);
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
