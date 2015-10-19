#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class DirectIntegratorEMS : public Integrator {
public:
    DirectIntegratorEMS(const PropertyList &props) {
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
            Vector3f d = ray.d;
            return distantsDisk->sampleL(d);
        }

        //get the Normal at the intersection point
        Vector3f n = Vector3f(its.shFrame.n);

        //get the Number of lights from the scene
        const  std::vector<Emitter *> lights = scene->getEmitters();
        uint32_t nLights = lights.size();


        //get the radiance of hitten object

        //cout << "meshEmitter = " << meshEmitter << endl;
        Color3f Le(0.0f, 0.0f, 0.0f);
        if (its.mesh->isEmitter() ) {
            const Emitter* areaLightEM = its.mesh->getEmitter();
            const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
            Le = aEM->sampleL(-ray.d, its.shFrame.n, its);
        }
        //cout << "Le = " << Le.toString() << endl;

            //get the asigned BSDF
            const BSDF* curBSDF = its.mesh->getBSDF();

            Color3f totalLight(0.0f, 0.0f, 0.0f);

            //iterated over all lights
            for (uint32_t var = 0; var < nLights; ++var) {

                //check if the light is a point lightsource
                int N = 1;
                Color3f totalLightRadiance(0.0f, 0.0f, 0.0f);

                //for every sample
                for (int i = 0; i < N; ++i) {
                    VisibilityTester vis;
                    Vector3f wi;
                    float pdf;
                    //create a query for the brdf
                    Vector3f wo = - ray.d;
                    Color3f Ld(0.0f, 0.0f, 0.0f);
                    Color3f f(0.0f, 0.0f, 0.0f);
                    const Point2f sample = sampler->next2D();
                    //sample the light
                    Ld = lights[var]->sampleL(its.p, Epsilon, sample, &wi, &pdf, &vis);
                    //std::cout << "Ld = " << Ld.toString() << std::endl;
                    //check if the light is visible
                    if(vis.Unoccluded(scene)) {

                        //calculate the abs value of n dot wi
                        float costerm = n.dot(wi);

                        //transform to the local frame
                        wo = its.toLocal(wo);
                        wi = its.toLocal(wi);

                        const BSDFQueryRecord query = BSDFQueryRecord(wi, wo, EMeasure::ESolidAngle);

                        //calculate the BSDF
                        f = curBSDF->eval(query);

                        //add the light up
                        if(pdf > 0.0)
                            totalLightRadiance += (Ld * f * std::abs(costerm)) / pdf;
                        //totalLight += (Ld * std::abs(costerm));

                    }
                }
                totalLightRadiance /= float(N);
                totalLight += totalLightRadiance;

            }
            return Le + totalLight;



    }

    std::string toString() const {
        return "DirectIntegratorEMS[]";
    }
private:

};

NORI_REGISTER_CLASS(DirectIntegratorEMS, "direct_ems");
NORI_NAMESPACE_END
