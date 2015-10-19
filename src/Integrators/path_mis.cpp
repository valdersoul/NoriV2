#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class PathIntegratorMIS : public Integrator {
public:
    PathIntegratorMIS (const PropertyList &props) {
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

        //get the Number of lights from the scene
        const  std::vector<Emitter *> lights = scene->getEmitters();
        uint32_t nLights = lights.size();

        Color3f tp(1.0f, 1.0f, 1.0f);
        Color3f L(0.0f, 0.0f, 0.0f);
        Ray3f pathRay(ray.o, ray.d);

        bool deltaFlag = true;

        while(true) {

            if (its.mesh->isEmitter() && deltaFlag) {
                const Emitter* areaLightEM = its.mesh->getEmitter();
                const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
                L += tp * aEM->sampleL(-pathRay.d, its.shFrame.n, its);
            }

            //Light sampling
            //randomly select a lightsource
            uint32_t var = uint32_t(std::min(sampler->next1D()*nLights, float(nLights) - 1.0f));

            //init the light color
            Color3f Li(0.0f, 0.0f, 0.0f);
            Color3f Ld(1.0f, 1.0f, 1.0f);

            //create a sample for the light
            const BSDF* curBSDF = its.mesh->getBSDF();
            const Point2f lightSample = sampler->next2D();
            VisibilityTester vis;
            Vector3f wo;
            float lightpdf;
            float bsdfpdf;
            Normal3f n = its.shFrame.n;

            deltaFlag = curBSDF->isDeltaBSDF();

            //sample the light

            {
                Li = lights[var]->sampleL(its.p, Epsilon, lightSample , &wo, &lightpdf, &vis);
                lightpdf /= float(nLights);
                //check if the pdf of the sample is greater than 0 and if the color is not black
                if(lightpdf > 0 && Li.maxCoeff() != 0.0f) {
                    //calculate the cosine term wi in my case the vector to the light
                    float cosTerm = std::abs(n.dot(wo));
                    const BSDFQueryRecord queryEM = BSDFQueryRecord(its.toLocal(- pathRay.d), its.toLocal(wo), EMeasure::ESolidAngle);
                    Color3f f = curBSDF->eval(queryEM);

                    if(f.maxCoeff() != 0.0f && vis.Unoccluded(scene)) {
                        bsdfpdf = curBSDF->pdf(queryEM);
                        float weight = BalanceHeuristic(float(1), lightpdf, float(1), bsdfpdf);
                        if(curBSDF->isDeltaBSDF())  weight = 1.0f;
                        if(bsdfpdf > 0.0f) {
                            Ld = (weight * f * Li * cosTerm) / lightpdf;
                            L += tp * Ld;
                        } else {
                            cout << "bsdfpdf = " << bsdfpdf  << endl;
                        }
                    }
                }
            }

            //Material part
            BSDFQueryRecord queryMats = BSDFQueryRecord(its.toLocal(-pathRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);

            Color3f fi =  curBSDF->sample(queryMats, sampler->next2D());
            bsdfpdf = curBSDF->pdf(queryMats);
            lightpdf = 0.0f;
            if(fi.maxCoeff() > 0.0f) {
                if(bsdfpdf > 0.0f) {
                    Ray3f shadowRay(its.p, its.toWorld(queryMats.wo));
                    Intersection lightIsect;

                     if (scene->rayIntersect(shadowRay, lightIsect)) {
                         if(lightIsect.mesh->isEmitter()){
                            const Emitter* areaLightEMcur = lightIsect.mesh->getEmitter();
                            const areaLight* aEMcur = static_cast<const areaLight *> (areaLightEMcur);

                            Li = aEMcur->sampleL(-shadowRay.d, lightIsect.shFrame.n, lightIsect);
                            lightpdf = aEMcur->pdf(its.p, (lightIsect.p - its.p).normalized(), lightIsect.p, Normal3f(lightIsect.shFrame.n));
                         }
                     } else {
                         const Emitter* distantsDisk = scene->getDistantEmitter();
                         if(distantsDisk != nullptr ) {
                             //check if THIS is right!
                             Li = distantsDisk->sampleL(lightIsect.toWorld(queryMats.wo));
                             lightpdf = distantsDisk->pdf(Point3f(0.0f), wo, Point3f(0.0f), Normal3f(0.0f));
                         }
                     }
                     lightpdf /= float(nLights);
                     //calculate the weights
                     float weight = BalanceHeuristic(float(1), bsdfpdf, float(1), lightpdf);


                     //check if the lightcolor is not black
                     if(Li.maxCoeff() > 0.0f  && lightpdf > 0.0f ) {
                         //wo in my case the vector to the light
                         Ld = weight * Li * fi;
                         L += tp * Ld;
                     }
                }



                 tp *= fi;
            } else {
                break;
            }

            wo = its.toWorld(queryMats.wo);
            pathRay = Ray3f(its.p, wo);

            if (!scene->rayIntersect(pathRay, its)) {
                const Emitter* distantsDisk = scene->getDistantEmitter();
                if(distantsDisk != nullptr ) {
                    //sample the distant disk light
                    Vector3f d = pathRay.d;
                    L += tp * distantsDisk->sampleL(d);
                }


                break;
            }

            float maxCoeff = tp.maxCoeff();
            float q = std::min(0.99f, maxCoeff);
            if(q < sampler->next1D()){
                break;
            }
            tp /= q;
        }

        return L;

    }

    static inline float BalanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        return (nf * fPdf) / (nf * fPdf + ng * gPdf);
    }

    std::string toString() const {
        return "PathIntegratorMIS []";
    }
private:

};

NORI_REGISTER_CLASS(PathIntegratorMIS , "path_mis");
NORI_NAMESPACE_END
