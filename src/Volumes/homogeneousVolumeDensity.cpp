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

#include <nori/warp.h>
#include <nori/volume.h>
#include <nori/bbox.h>


NORI_NAMESPACE_BEGIN


class HomogeneousVolumeDensity : public Volume {
public:
    HomogeneousVolumeDensity(const PropertyList &propList) :
        Volume(propList.getTransform("toWorld", Transform()))
    {
        sig_a = propList.getColor("sig_a", Color3f(0.5f));
        sig_s = propList.getColor("sig_s", Color3f(0.5f));
        sig_t = sig_a + sig_s;
        albido = sig_s / sig_t;
        le = propList.getColor("le", Color3f(0.5f));
        Vector3f BBMin = propList.getVector("BBMin", Vector3f(0.f, 0.f, 0.f));
        Vector3f BBMax = propList.getVector("BBMax", Vector3f(0.f, 0.f, 0.f));
        extent = BoundingBox3f(BBMin, BBMax);

    }

    bool IntersectP(const Ray3f &ray, float &t1, float &t2) const {
        Ray3f localRay = WorldToVolume * ray;
        return extent.rayIntersect(localRay, t1, t2);
    }
    /*
    BoundingBox3f WorldBound() const {
        return VolumeToWorld
    }
    */

    Color3f sigma_a(const Point3f &p, const Vector3f &v) const {
        return extent.contains(WorldToVolume * p) ? sig_a : Color3f(0.f);
    }

    Color3f sigma_s(const Point3f &p, const Vector3f &v) const {
        return sig_s;
    }

    Color3f sigma_t(const Point3f &p, const Vector3f &v) const {
        if (extent.contains(WorldToVolume * p)) {
            return (sig_a + sig_s);
        }
        return 0.;
    }

    Color3f getSimga_t () const {
        return sig_t;
    }

    Color3f getSimga_s () const {
        return sig_s;
    }

    Color3f Lve(const Point3f &p, const Vector3f &v) const {
        return extent.contains(WorldToVolume * p) ? le : Color3f(0.f);
    }
    //NOT USED
    float pdf(const Point3f &p, const Vector3f &wi, const Vector3f &wo) const {
        PhaseFunctionQueryRecord phaseQry(wi, wo, ESolidAngle);

        return m_phase->eval(phaseQry);
    }


    Color3f tau(const Ray3f &ray, float step = 1.f, float offset = 0.5) const {
        float t1;
        float t2;
        if(!IntersectP(ray, t1, t2)) return Color3f(0.f);
        return (t2 - t1) * (sig_a + sig_s);
    }

    void samplePhaseFunction(PhaseFunctionQueryRecord &pQry, const Point2f &sample) const {
        m_phase->sample(pQry, sample);
    }

    float getFreeFlightDistance(float sample, int channel, float &pdfMedium, float &pdfSurface) const  {       
        Color3f sigT = getSimga_t();
        //sample /= 2.0f;
        float d = - logf(1.0f - sample) / sigT(channel);

        return d;

    }
    float getTransmittance(float distance, int channel) const {
        Color3f sigT = getSimga_t();
        float T = exp(-sigT(channel) * distance);
        return T;
    }

    Color3f getTransmittance(float distance) const {
        Color3f temp = (-distance) * getSimga_t();
        return temp.expElemwise();
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "HomogeneousVolumeDensity[\n"
            "  sig_a = %s\n"
            "  sig_s = %s\n"
            "  le = %s\n"
            "]", sig_a.toString(), sig_s.toString(), le.toString());
    }


    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EPhaseFunction:
                if (m_phase)
                    throw NoriException(
                        "Mesh: tried to register multiple Phase functions!");
                m_phase = static_cast<PhaseFunction *>(obj);
                break;

            default:
                throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }



    PhaseFunction* getPhaseFkt() const {
        return m_phase;
    }

    float getTransmittancePDF(float distance, int channel) const {
        Color3f sigT = getSimga_t();
        float pdf = sigT[channel] * exp(-sigT[channel] * distance);
        return pdf;
    }

private:
    PhaseFunction *m_phase = nullptr;
    Color3f sig_a;
    Color3f sig_s;
    Color3f sig_t;
    Color3f albido;
    Color3f le;
    BoundingBox3f extent;
};

NORI_REGISTER_CLASS(HomogeneousVolumeDensity, "homogeneousVolumeDensity");
NORI_NAMESPACE_END
