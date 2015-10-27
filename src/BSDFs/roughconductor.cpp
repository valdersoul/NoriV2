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
#include <nori/metalConsts.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/// Ideal RoughConductor BRDF
class RoughConductor : public BSDF {
public:
    RoughConductor(const PropertyList &propList) {
        materialName = propList.getString("materialName", "");

        if(materialName != "") {
            // load the specific material properties
            if(!getMaterialProperties(materialName, m_eta, m_k)){
                std::cout << "Material " << materialName << " not found!" << std::endl;
            }
        } else {
            m_eta = propList.getColor("eta");
            m_k = propList.getColor("k");
        }

        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);
    }

    Color3f fresnelCond(const float cosTheta) const {
        // compute the fresnel coefficients for all channels
        return Color3f(
                conductorReflectance(m_eta.x(), m_k.x(), cosTheta),
                conductorReflectance(m_eta.y(), m_k.y(), cosTheta),
                conductorReflectance(m_eta.z(), m_k.z(), cosTheta)
            );
    }


    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &bRec) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &) const {
        if (Frame::cosTheta(bRec.wi) <= 0) 
            return Color3f(0.0f);

        return Color3f(0.0f);
    }

    std::string toString() const {
        return tfm::format(
                    "RoughConductor[\n"
                    "  material = %s,\n"
                    "  m_eta = %f,\n"
                    "  m_k = %f\n",
                    "  m_alpha = %f,\n"
                    "]",
                    materialName, m_eta, m_k, m_alpha);
    }
    bool isDeltaBSDF() const {
        return true;
    }
private:
    std::string materialName;
    //index of refraction
    Color3f m_eta;

    // index of absorbtion
    Color3f m_k;

    // the roughness of the material
    float m_alpha;
};

NORI_REGISTER_CLASS(RoughConductor, "roughConductor");
NORI_NAMESPACE_END
