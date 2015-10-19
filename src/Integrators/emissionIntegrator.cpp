#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>
#include <nori/volume.h>
#include <nori/color.h>
#include <math.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class EmissionIntegrator : public VolumeIntegrator {
public:
    EmissionIntegrator(const PropertyList &props) {
        stepSize = props.getFloat("stepSize", 0.5f);
    }
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        return Color3f(0.0f);
    }
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, Color3f &T) const {
        const std::vector<Volume *> volumes = scene->getVolumes();

        Color3f result(0.0f, 0.0f, 0.0f);
        if(volumes.size() == 0) {
            T = Color3f(1.0f, 1.0f, 1.0f);
            return result;
        }

        for (unsigned int i = 0; i < volumes.size(); ++i) {
            Volume *vr = volumes[i];
            float t0 = 0.0f ;
            float t1 = 0.0f ;
            if(!vr || !vr->IntersectP(ray, t0, t1) || (t1 -t0) == 0.0f){
                //CHECK WHAT TODO
                //T = Color3f(1.0f, 1.0f, 1.0f);

                continue;
            }

            Color3f Lv(0.0f, 0.0f, 0.0f);

            int nSamples = int(ceil((t1 - t0) / stepSize));
            float step = (t1 - t0) / nSamples;
            Color3f Tr(1.0f, 1.0f, 1.0f);
            Point3f p = ray.o;
            Point3f pPrev;
            Vector3f w = -ray.d;
            t0 += sampler->next1D() * step;
            for (int i = 0; i < nSamples; ++i, t0 += step) {
                // Advance to sample at _t0_ and update _T_
                pPrev = p;
                p = ray.o + (ray.d * t0);
                Ray3f tauRay(pPrev, p - pPrev, 0.0f, ray.maxt);
                Color3f stepTau = vr->tau(tauRay, .5f * stepSize, sampler->next1D());
                stepTau = - stepTau;
                Color3f curTr = stepTau.expElemwise();
                Tr = Tr * curTr;

                // Possibly terminate ray marching if transmittance is small
                if (Tr.y() < 1e-3) {
                    const float continueProb = .5f;
                    if (sampler->next1D() > continueProb) {
                        Tr = Color3f(0.0f, 0.0f, 0.0f);
                        break;
                    }
                    Tr /= continueProb;
                }

                // Compute emission-only source term at _p_
                Lv += Tr * vr->Lve(p, w);
            }
            T = Tr;
            result += Lv * step;
        }

        return result;
    }

    Color3f Transmittance(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        const std::vector<Volume *> volumes = scene->getVolumes();
        if(volumes.size() == 0)
            return Color3f(1.0f, 1.0f, 1.0f);

        Color3f tau(0.0f, 0.0f, 0.0f);

        for (unsigned int i = 0; i < volumes.size(); ++i) {
            float offset = sampler->next1D();
            tau += volumes[i]->tau(ray, stepSize, offset);
        }
        tau = - tau;
        Color3f trans = tau.expElemwise();
        return trans;
    }

    std::string toString() const {
        return "EmissionIntegrator[]";
    }
private:
    float stepSize;
};

NORI_REGISTER_CLASS(EmissionIntegrator, "emissionIntegrator");
NORI_NAMESPACE_END
