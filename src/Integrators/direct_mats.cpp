#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class DirectIntegratorMATS : public Integrator {
public:
    DirectIntegratorMATS(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;

        //check if the ray intersects the scene
        if (!scene->rayIntersect(ray, its)) {
            //check if a distant disk light is set
            const Emitter* distantsDisk = scene->getDistantEmitter();
            if(distantsDisk == nullptr ) return Color3f(0.0f);

            //sample the distant disk light
            return distantsDisk->sampleL(ray.d);
        }

        //get the radiance of hitten object
        Color3f Le(0.0f, 0.0f, 0.0f);
        if (its.mesh->isEmitter()  ) {
            const Emitter* areaLightEM = its.mesh->getEmitter();
            const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
            Le = aEM->sampleL(-ray.d, its.shFrame.n, its);
        }

        //get the asigned BSDF
        const BSDF* curBSDF = its.mesh->getBSDF();

        Color3f Ld(0.0f, 0.0f, 0.0f);
        Color3f f(0.0f, 0.0f, 0.0f);
        Color3f totalLight(0.0f, 0.0f, 0.0f);

        //transform to the local frame
        //create a BRDF Query
        BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-ray.d), Vector3f(0.0f), EMeasure::ESolidAngle);

        //sample the BRDF
        Color3f mats =  curBSDF->sample(query, sampler->next2D());

        if(mats.maxCoeff() > 0.0f) {
            //Check for the light source
            Vector3f wo = its.toWorld(query.wo);
            Ray3f shadowRay(its.p, wo);
            Intersection itsShadow;
            if (scene->rayIntersect(shadowRay, itsShadow)) {
                //intersection check if mesh is emitter
                if(itsShadow.mesh->isEmitter()){
                    Ld = itsShadow.mesh->getEmitter()->radiance();
                }
            } else {
                //check for distant disk light
                const Emitter* distantsDisk = scene->getDistantEmitter();
                if(distantsDisk != nullptr ) Ld = distantsDisk->sampleL(wo);
            }

            totalLight += Ld * mats;
        }

        return Le + totalLight;
    }

    std::string toString() const {
        return "DirectIntegratorMATS[]";
    }
private:

};

NORI_REGISTER_CLASS(DirectIntegratorMATS, "direct_mats");
NORI_NAMESPACE_END
