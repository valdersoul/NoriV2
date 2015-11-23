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

#include <nori/medium.h>

NORI_NAMESPACE_BEGIN


class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium   (const PropertyList &propList)
    {

    }

    bool sampleDistance(const Ray3f &ray, MediumSamplingRecord &mRec, Sampler *sampler) const {
        return false;
    }

    /**
     * \brief Compute the 1D density of sampling distance \a ray.maxt
     * along the ray using the sampling strategy implemented by
     * \a sampleDistance.
     *
     * The function computes the continuous densities in the case of
     * a successful \ref sampleDistance() invocation (in both directions),
     * as well as the Dirac delta density associated with a failure.
     * For convenience, it also stores the transmittance along the
     * supplied ray segment within \a mRec.
     */
    void eval(const Ray3f &ray, MediumSamplingRecord &mRec) const {

    }
    /**
     * \brief Compute the transmittance along a ray segment
     *
     * Computes the transmittance along a ray segment
     * [mint, maxt] associated with the ray. It is assumed
     * that the ray has a normalized direction value.
     */
    Color3f evalTransmittance(const Ray3f &ray,
        Sampler *sampler = NULL) const {
        return Color3f(0.0f);
    }

    /// Determine whether the medium is homogeneous
    bool isHomogeneous() const  {
        return true;
    }

    /// Add a child
    void addChild(NoriObject *obj) {}

    /// Return a string representation
    std::string toString() const {
        return "";
    }

};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogeneousMedium");
NORI_NAMESPACE_END
