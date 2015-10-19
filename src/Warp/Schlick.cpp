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

#include <nori/phaseFunction.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Diffuse / Lambertian BRDF model
 */
class Schlick : public PhaseFunction {
public:
    Schlick(const PropertyList &propList) {
        m_k = propList.getFloat("k", 0.5f);
        //m_k = 1.55f * g - 0.55 * g * g * g;
    }


    /// Compute the density of \ref sample() wrt. solid angles
    float eval(PhaseFunctionQueryRecord &bRec) const {
        bRec.pdf  = Warp::squareToSchlickPDF(bRec.wo, m_k); //Warp::squareToSchlickPDF(bRec.wo, bRec.wo, m_k);
        return 1.0f;
    }

    /// Draw a a sample from the BRDF model
    float sample(PhaseFunctionQueryRecord &bRec, const Point2f &sample) const {
        bRec.wo = Warp::squareToSchlick(sample, m_k);
        eval(bRec);

        return 1.0f;
    }


    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Schlick[\n"
            "  k = %f\n"
            "]", m_k);
    }


private:
    float m_k;
};

NORI_REGISTER_CLASS(Schlick, "schlick");
NORI_NAMESPACE_END
