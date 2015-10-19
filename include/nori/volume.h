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

#if !defined(__NORI_VOLUME_H)
#define __NORI_VOLUME_H

#include <nori/object.h>
#include <nori/integrator.h>
#include <nori/phaseFunction.h>

NORI_NAMESPACE_BEGIN



/**
 * \brief Superclass of all bidirectional scattering distribution functions
 */
class Volume : public NoriObject {
public:
    Volume(const Transform &l2v) :
        VolumeToWorld(l2v),
        WorldToVolume(l2v.inverse())
    {
        //WARN if matrix contains scale
    }
    Volume() :
        VolumeToWorld(Transform()),
        WorldToVolume(Transform())
    {

    }
    virtual bool IntersectP(const Ray3f &ray, float &t1, float &t2) const = 0;

    //virtual BoundingBox3f WorldBound() const = 0;
    virtual Color3f getSimga_t() const = 0;

    virtual Color3f getSimga_s() const = 0;

    virtual Color3f sigma_a(const Point3f &p, const Vector3f &v) const = 0;

    virtual Color3f sigma_s(const Point3f &p, const Vector3f &v) const = 0;

    virtual Color3f sigma_t(const Point3f &p, const Vector3f &v) const = 0;

    virtual Color3f Lve(const Point3f &p, const Vector3f &v) const = 0;

    virtual float pdf(const Point3f &p, const Vector3f &wi, const Vector3f &wo) const = 0;

    virtual Color3f tau(const Ray3f &ray, float step = 1.f, float offset = 0.5) const = 0;

    virtual void samplePhaseFunction(PhaseFunctionQueryRecord &pQry, const Point2f &sample) const = 0;

    virtual float getFreeFlightDistance(float sample, int channel, float &pdfMedium, float &pdfSurface) const = 0;
    virtual float getTransmittance(float distance, int channel) const = 0;
    virtual Color3f getTransmittance(float distance) const = 0;
    virtual float getTransmittancePDF(float distance, int channel) const = 0;
    virtual PhaseFunction* getPhaseFkt() const = 0;

    virtual void addChild(NoriObject *obj) {}

    EClassType getClassType() const { return EMedium; }
protected:
     //Variabeles
    const Transform VolumeToWorld;
    const Transform WorldToVolume;

};

class VolumeIntegrator : public Integrator {
public:
    // VolumeIntegrator Interface
    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, Color3f &T) const = 0;
    virtual Color3f Transmittance(const Scene *scene, Sampler *sampler, const Ray3f &ray) const = 0;
};

NORI_NAMESPACE_END

#endif
