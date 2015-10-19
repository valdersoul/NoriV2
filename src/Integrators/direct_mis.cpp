#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class DirectIntegratorMIS : public Integrator {
public:
    DirectIntegratorMIS(const PropertyList &props) {
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
        Normal3f n = its.shFrame.n;

        //get the Number of lights from the scene
        const  std::vector<Emitter *> lights = scene->getEmitters();
        uint32_t nLights = lights.size();

        //init the light terms
        Color3f Le(0.0f, 0.0f, 0.0f);
        Color3f Ld(0.0f, 0.0f, 0.0f);

        //check if the mesh is an emitter
        if (its.mesh->isEmitter())  {
            const Emitter* areaLightEM = its.mesh->getEmitter();
            const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
            Le = aEM->sampleL(-ray.d, its.shFrame.n, its);
        }

        //get the asigned BSDF
        const BSDF* curBSDF = its.mesh->getBSDF();

        uint32_t var = uint32_t(sampler->next1D()*nLights);
        //iterated over all lights
        //for (u_int32_t var = 0; var < nLights; ++var)
        {

            //init the light color
            Color3f Li(0.0f, 0.0f, 0.0f);

            //create a sample for the light
            const Point2f lightSample = sampler->next2D();
            VisibilityTester vis;
            Vector3f wo;
            Vector3f wi = - ray.d;
            float lightpdf;
            float bsdfpdf;

            //sample the light
            Li = lights[var]->sampleL(its.p, Epsilon, lightSample , &wo, &lightpdf, &vis);
            lightpdf /= float(nLights);
            //check if the pdf of the sample is greater than 0 and if the color is not black
            if(lightpdf > 0 && Li.maxCoeff() != 0.0f) {

                //calculate the cosine term wi in my case the vector to the light
                float cosTerm = std::abs(n.dot(wo));

                //calculate the BSDF
                const BSDFQueryRecord queryEM = BSDFQueryRecord(its.toLocal(wi), its.toLocal(wo), EMeasure::ESolidAngle);
                Color3f f = curBSDF->eval(queryEM);

                //check if the color is not black and if the scene is unoccluded
                if(f.maxCoeff() != 0.0f && vis.Unoccluded(scene)) {

                    //check if the light is a delta light
                    if(lights[var]->isDeltaLight() ){
                        Ld += (f * Li * cosTerm) / lightpdf;
                    } else {
                        //get the pdf of the BSDF and compute the weights
                        bsdfpdf = curBSDF->pdf(queryEM);
                        float weight = BalanceHeuristic(float(1), lightpdf, float(1), bsdfpdf);
                        if(bsdfpdf > 0.0f) {
                            Ld += (weight * f * Li * cosTerm) / lightpdf;
                        } else {
                            cout << "bsdfpdf = " << bsdfpdf  << endl;
                        }
                    }
                }
            }

            //sample the BSDF only if the the light is no delta light -> pdf allways zero
            if(!lights[var]->isDeltaLight()){

                //create a BSDF query, sample, eval it and get the pdf of the sample
                BSDFQueryRecord queryMats = BSDFQueryRecord(its.toLocal(-ray.d), Vector3f(0.0f), EMeasure::ESolidAngle);
                Color3f mats =  curBSDF->sample(queryMats, sampler->next2D());
                bsdfpdf = curBSDF->pdf(queryMats);

                //check if the sampled direction intersects with a lightsource
                lightpdf = 0.0f;
                Intersection lightIsect;
                Color3f Li(0.0f);
                Ray3f shadowRay(its.p, its.toWorld(queryMats.wo));

                //check if the pdf of the sample is greater than 0 and if the color is not black
                if(bsdfpdf > 0.0f && mats.maxCoeff() > 0.0f) {

                    //check for intersections
                    if (scene->rayIntersect(shadowRay, lightIsect)) {
                        //intersection check if mesh is emitter
                        if(lightIsect.mesh->isEmitter()){
                            Li = lightIsect.mesh->getEmitter()->radiance();
                            lightpdf = lightIsect.mesh->getEmitter()->pdf(its.p, (lightIsect.p - its.p).normalized(), lightIsect.p, Normal3f(lightIsect.shFrame.n));
                        }
                    } else {
                        //check for distant disk light
                        const Emitter* distantsDisk = scene->getDistantEmitter();
                        if(distantsDisk != nullptr ) {
                            //check if THIS is right!
                            Li = distantsDisk->sampleL(its.toWorld(queryMats.wo));
                            lightpdf = distantsDisk->pdf(Point3f(0.0f), wo, Point3f(0.0f), Normal3f(0.0f));                            
                        }
                    }
                    lightpdf /= float(nLights);
                    //calculate the weights
                    float weight = BalanceHeuristic(float(1), bsdfpdf, float(1), lightpdf);


                    //check if the lightcolor is not black
                    if(Li.maxCoeff() > 0.0f  && lightpdf > 0.0f ) {
                        //wo in my case the vector to the light

                        Ld += weight * Li * mats;
                    }
                }
            }
        } //end for

        return Le + Ld /* float(nLights)*/ ;

    }

    static inline float BalanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        return (nf * fPdf) / (nf * fPdf + ng * gPdf);
    }

    std::string toString() const {
        return "DirectIntegratorMIS[]";
    }
private:

};

NORI_REGISTER_CLASS(DirectIntegratorMIS, "direct_mis");
NORI_NAMESPACE_END
