#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>


#include <iostream>
#include <stdint.h>

NORI_NAMESPACE_BEGIN

class testtexture : public Integrator {
public:
    testtexture(const PropertyList &props) {

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
        Normal3f n = Vector3f(its.shFrame.n);

        //check if the mesh has a bumpmap

        //get the Number of lights from the scene
        const  std::vector<Emitter *> lights = scene->getEmitters();
        uint32_t nLights = lights.size();


        //get the radiance of hitten object
        const Emitter* meshEmitter = its.mesh->getEmitter();
        //cout << "meshEmitter = " << meshEmitter << endl;
        Color3f Le(0.0f, 0.0f, 0.0f);
        if (meshEmitter != nullptr  ) Le = meshEmitter->radiance();
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

                        //check if the mesh has texture
                        Color3f tex(0.0f, 0.0f, 0.0f);
                        if(its.mesh->hasTexture()){
                            if((!isnan(its.uv(0))) && (!isnan(its.uv(1))))
                                tex = its.mesh->getTexture()->sample(its.uv);
                        } else {
                            tex = f;
                        }

                        //add the light up
                        if(pdf > 0.0)
                            totalLightRadiance += (Ld * tex * std::abs(costerm)) / pdf;
                            //totalLightRadiance += (Ld * f * std::abs(costerm)) / pdf;


                    }
                }
                totalLightRadiance /= float(N);
                totalLight += totalLightRadiance;

            }
            return Le + totalLight;
    }

    std::string toString() const {
        return "testtexture []";
    }
private:

};

NORI_REGISTER_CLASS(testtexture, "testtexture");
NORI_NAMESPACE_END
