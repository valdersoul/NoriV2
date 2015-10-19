#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectIntegrator : public Integrator {
public:
    DirectIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        //get the Normal at the intersection point
        Vector3f n = Vector3f(its.shFrame.n);

        //get the Number of lights from the scene

        const  std::vector<Emitter *> lights = scene->getEmitters();
        uint32_t nLights = lights.size();


        //check for lights
        if(nLights > 0) {

            //get the asigned BSDF
            const BSDF* curBSDF = its.mesh->getBSDF();
            if(curBSDF) {
                Color3f totalLight(0.0f);
                //uterated over all lights
                //std::cout << "start light comp" << std::endl;
                for (uint32_t var = 0; var < nLights; ++var) {

                    VisibilityTester vis;
                    Vector3f wi;
                    float pdf;
                    const Point2f ls;

                    //sample the light
                    Color3f Ld =lights[var]->sampleL(its.p, Epsilon, ls, &wi, &pdf, &vis);

                    //check if the light is visible
                    if(vis.Unoccluded(scene)) {

                        //create a query for the brdf
                        Vector3f wo = - ray.d;

                        //calculate the abs value of n dot wi
                        float costerm = n.dot(wi);

                        //transform to the local frame
                        wo = its.toLocal(wo);
                        wi = its.toLocal(wi);

                        const BSDFQueryRecord query = BSDFQueryRecord(wi, wo, EMeasure::ESolidAngle);

                        //calculate the BSDF
                        Color3f f = curBSDF->eval(query);

                        //add the light up
                        totalLight += (Ld * f * std::abs(costerm));

                    }
                }

                return totalLight;
            }
        }
        return Color3f(0.0f);


    }

    std::string toString() const {
        return "DirectIntegrator[]";
    }
private:

};

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END
