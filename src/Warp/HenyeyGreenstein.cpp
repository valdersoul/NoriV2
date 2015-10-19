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
class HenyeyGreenstein : public PhaseFunction {
public:
    HenyeyGreenstein(const PropertyList &propList) {
        m_g = propList.getFloat("g", 0.5f);
    }


    /// Compute the density of \ref sample() wrt. solid angles
    float eval(PhaseFunctionQueryRecord &bRec) const {
        bRec.pdf  = Warp::squareToHenyeyGreensteinPDF(bRec.wo, m_g);
        return 1.0f;
    }

    /// Draw a a sample from the BRDF model
    float sample(PhaseFunctionQueryRecord &bRec, const Point2f &sample) const {
        bRec.wo = Warp::squareToHenyeyGreenstein(sample, m_g);
        eval(bRec);

        return 1.0f;
    }


    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "HenyeyGreenstein[\n"
            "  g = %f\n"
            "]", m_g);
    }


private:
    float m_g;
};

NORI_REGISTER_CLASS(HenyeyGreenstein, "henyeyGreenstein");
NORI_NAMESPACE_END
