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

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {

        bRec.measure = EDiscrete;

        //check if your are inside or not
        //compute the angle between the normal and the incoming direction
        float theta1 = Frame::cosTheta(bRec.wi);

        float fresnelTerm = fresnel(theta1, m_extIOR, m_intIOR);

        //check if entering or leaving
        bRec.entering = (theta1 > 0);


        float et = m_extIOR;
        float ei = m_intIOR;
        Vector3f n = Vector3f(0.0, 0.0, 1.0);
        if(!bRec.entering) {
             std::swap(ei, et);
             n(2) = -1.0;
        }


        // check you should a reflection or refraction
        if(sample(0) < fresnelTerm) {
            //mirror like reflection
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                 bRec.wi.z()
            );
            bRec.eta = 1.0f;
        } else {
            float snellsTerm = et / ei;

            float wdotn = bRec.wi.dot(n);

            Vector3f fstPart = (bRec.wi - (wdotn) * n);
            float srqtTerm = std::sqrt(1.0f - snellsTerm * snellsTerm * (1.0f - wdotn * wdotn));

            bRec.wo = - snellsTerm * fstPart - n * srqtTerm;
            bRec.eta = bRec.entering?  m_intIOR / m_extIOR : m_extIOR / m_intIOR;
        }

        return Color3f(1.0f);
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }

    bool isDeltaBSDF() const {
        return true;
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
