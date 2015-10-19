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

#if !defined(__NORI_EMITTER_H)
#define __NORI_EMITTER_H

#include <nori/object.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/photon.h>
#include <nori/bitmap.h>
#include <filesystem/resolver.h>

NORI_NAMESPACE_BEGIN

/**
 * @brief The VisibilityTester struct
 * tests if a there is any object between to points
 */
struct VisibilityTester {
    Ray3f r;
    /**
     * @brief SetSegment the the ray
     * @param p1 the first point of the ray (origin)
     * @param eps1 a offset value for the first point
     * @param p2 the second point of the ray
     * @param eps2 a offset value for the second point
     */
    void SetSegment(const Point3f &p1, float eps1,
                    const Point3f &p2, float eps2 /*, float time*/) {
        //calculate the direction
        Vector3f d = (p2 - p1);

        //calculate the distance
        float dist = (p2 - p1).norm();

        //create the ray segment
        r = Ray3f(p1, d / dist, eps1 , dist * (1.f - eps2)); /*, time);*/
    }

    /**
     * @brief SetRay create a ray for the VisibilityTester
     * @param p the origin
     * @param w the direction (should be normalized)
     */
    void SetRay(const Point3f &p, const Vector3f &w/* , float eps, float time*/) {
        r = Ray3f(p, w);
    }
    /**
     * @brief Unoccluded check if the origin is visible
     * @param scene the current scene
     * @return true if the the point is not occluded
     */
    bool Unoccluded(const Scene *scene) const {
        return !scene->rayIntersect(r);
    }
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:
    //Variables
    //the number of the samples used to sample the light
    int nSamples;

    /**
     * @brief Emitter Constructor
     * @param l2w Transformation from light to world space
     * @param nrSmaples the number of samples used to sample the light source
     */
    Emitter(const Transform &l2w, const int nrSmaples = 1) :
        nSamples(std::max(1, nrSmaples)),
        LightToWorld(l2w),
        WorldToLight(l2w.inverse())
    {
        //WARN if matrix contains scale
    }

    /**
     * @brief Emitter Constructor
     */
    Emitter() :
        nSamples(1),
        LightToWorld(Transform()),
        WorldToLight(Transform())
    {

    }


    /// Release all memory
    virtual ~Emitter() { }
    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
    /**
     * @brief Power return the total power of the light
     * @return the total power
     */
    virtual Color3f Power(const Scene *) const = 0;

    virtual void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const = 0;


    /**
     * @brief sampleL sample the light source
     * @param p the current point of interessed
     * @param pEpsilon offset from the point p
     * @param wi direction to the light source
     * @param pdf the probalility of the sample
     * @param vis the initalized visiblity tester
     * @return a color sample of the light
     */
    virtual Color3f sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /* float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const = 0;

    /**
     * @brief sampleL sample the light source (getting Le)
     * @param p the current point of interessed
     */
    //virtual Color3f sampleL(Point3f &p) const = 0;

    //TEMPORY
    virtual Color3f sampleL(const Vector3f &d) const = 0;


    /**
     * @brief isDeltaLight check if the light is a delta light
     * @return true is the ligth is a delta light
     */
    virtual bool isDeltaLight() const = 0;

    virtual float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const = 0;

    void setMaxRadius(float radius){

    }

    virtual Color3f radiance() const = 0;
protected:
     //Variabeles
    const Transform LightToWorld;
    const Transform WorldToLight;


};

NORI_NAMESPACE_END

#endif /* __NORI_EMITTER_H */
