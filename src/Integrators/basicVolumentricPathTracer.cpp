#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>
#include <nori/volume.h>


#include <iostream>

NORI_NAMESPACE_BEGIN

class BasicPathIntegrator : public Integrator {
public:
    BasicPathIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }
    bool checkVolIntersetcion(Ray3f &pathRay, float &t0, float &t1, const std::vector<Volume *> &volumes, Volume *&curVol) const {
        curVol = nullptr;
        //TODO: check what happens if multiple volumens after each other
        for (unsigned int i = 0; i < volumes.size(); ++i) {
            //use the first one
            if(volumes[i]->IntersectP(pathRay, t0, t1)){
                curVol = volumes[i];
                if(t0 < 0) {
                    //inside the medium
                    t0 = 0.f;
                }

                return true;
            }
        }
        return false;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        bool debug = false;
        Color3f tp(1.0f, 1.0f, 1.0f);
        Color3f L(0.0f, 0.0f, 0.0f);
        Ray3f pathRay(ray.o, ray.d);
        const std::vector<Volume *> volumes = scene->getVolumes();
        float pdfMedium;
        float pdfSurface;

        while(true) {
            //check if there is a volume ray intersection
            Volume *curVol = nullptr;
            float t0 = 0.0f;
            float t1 = 0.0f;
            float ffDist = 0.0f;
            bool surIntersection = scene->rayIntersect(pathRay, its);

            if(checkVolIntersetcion(pathRay, t0, t1, volumes, curVol)) {
                int sigIndex = int(std::min(sampler->next1D() * 3.0f, 2.0f));
                ffDist = curVol->getFreeFlightDistance(sampler->next1D(), sigIndex, pdfMedium, pdfSurface);
            }


            //check if a volume exists, if the free fliegt distance is
            //inside the media and befor the next surface
            bool media = false;
            float volEnd = t1 - t0;
            if (curVol != nullptr &&   /* check if there is a volume && t0 <= ffDist */
                ffDist <= volEnd &&      /* check if the free flight distance is inside the medium */
                (!surIntersection || /* no mesh in that direction */
                 ffDist < (its.t - t0)))        /* intersection is after the volume */
            {
                media = true;
                //compute the postion
                Point3f o = pathRay(ffDist + t0);

                //Create a frame to transform to volume
                Frame w2v(pathRay.d);

                //create a query for the phase function
                PhaseFunctionQueryRecord pQry(w2v.toLocal(-pathRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);

                Color3f albido =  curVol->getSimga_s() / curVol->getSimga_t();
                if(albido.maxCoeff() > 0.0f) {
                    tp *= albido;
                } else {
                    break;
                }
                //sample the phase function
                curVol->samplePhaseFunction(pQry, sampler->next2D());
                Vector3f wo = w2v.toWorld(pQry.wo);
                pathRay = Ray3f(o, wo, 0.0, std::numeric_limits<float>::infinity());

            } else {
                //ray escapes the scene
                //check if the ray intersects the scene
                if (!surIntersection) {
                    //TODO: ADD MEDIA

                    //check if a distant disk light is set
                    const Emitter* distantsDisk = scene->getDistantEmitter();
                    if(distantsDisk == nullptr ) return Color3f(0.0f);

                    //sample the distant disk light
                    return tp * distantsDisk->sampleL(pathRay.d);
                }

                //get the radiance of hitten object
                if (its.mesh->isEmitter() ) {
                    const Emitter* areaLightEM = its.mesh->getEmitter();
                    const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
                    L += tp * aEM->sampleL(-pathRay.d, its.shFrame.n, its);
                }

                //get the asigned BSDF
                const BSDF* curBSDF = its.mesh->getBSDF();

                //transform to the local frame
                BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-pathRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);

                //sample the BRDF
                Color3f fi =  curBSDF->sample(query, sampler->next2D());

                //check for black brdf
                if(fi.maxCoeff() > 0.0f) {
                    tp *= fi;
                } else {
                    //stop
                    // hit a black brdf
                    break;
                }
                Vector3f wo = its.toWorld(query.wo);
                pathRay = Ray3f(its.p, wo);
            }

            //stop critirium russian roulette
            float maxCoeff = tp.maxCoeff();
            float q = std::min(0.99f, maxCoeff);
            if(q < sampler->next1D()){
                if(media && debug)
                    cout << "breaking"  << endl;
                break;
            }
            tp /= q;

        }
        return L;
    }

    std::string toString() const {
        return "BasicPathIntegrator[]";
    }

};

NORI_REGISTER_CLASS(BasicPathIntegrator, "basicPathIntegrator");
NORI_NAMESPACE_END
