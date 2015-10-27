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

NORI_NAMESPACE_BEGIN

/// Ideal Conductor BRDF
class Conductor : public BSDF {
public:
    Conductor(const PropertyList &propList) {
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
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        // check if the given vector pair is a perfect reflection
        Vector3f reflectionVec = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());

        // check if the reflection vector is the same as the incoming vector
        // take rounding errors into account
        if(std::abs(reflectionVec.dot(bRec.wo) - 1.0f) > DeltaEpsilon) {
            return Color3f(0.0f);
        }

        return fresnelCond(Frame::cosTheta(bRec.wi));
    }

    float pdf(const BSDFQueryRecord &bRec) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // check if the given vector pair is a perfect reflection
        Vector3f reflectionVec = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());

        // check if the reflection vector is the same as the incoming vector
        // take rounding errors into account
        if(std::abs(reflectionVec.dot(bRec.wo) - 1.0f) > DeltaEpsilon) {
            return 1.0f;
        }

        return 1.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &) const {
        if (Frame::cosTheta(bRec.wi) <= 0) 
            return Color3f(0.0f);

        // Reflection in local coordinates
        bRec.wo = Vector3f(
            -bRec.wi.x(),
            -bRec.wi.y(),
             bRec.wi.z()
        );
        bRec.measure = EDiscrete;

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        return fresnelCond(Frame::cosTheta(bRec.wi));
    }

    std::string toString() const {
        return tfm::format(
                    "conductor[\n"
                    "  material = %s,\n"
                    "  m_eta = %s,\n"
                    "  m_k = %s\n"
                    "]",
                    materialName, m_eta.toString(), m_k.toString());
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
};

NORI_REGISTER_CLASS(Conductor, "conductor");
NORI_NAMESPACE_END
