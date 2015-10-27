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

#if !defined(__NORI_MICROFACET_H)
#define __NORI_BSDF_H

#include <nori/object.h>
#include <nori/warp.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Class for a mircofacet
 */
class Microfacet {

public:
    static float D(float alpha, const Vector3f &m);
    static float G(float alpha,const Vector3f &i, const Vector3f &o, const Vector3f &m);
    static float pdf(float alpha, const Vector3f &m);
    static Vector3f sample(float alpha, Point2f xi);
private:
    static float G1(float alpha, const Vector3f &v, const Vector3f &m);
    static float CHIPlus(const float c);
};

NORI_NAMESPACE_END

#endif /* __NORI_MICROFACET_H */
