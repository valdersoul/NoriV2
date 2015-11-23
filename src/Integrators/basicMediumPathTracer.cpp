#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>
#include <nori/medium.h>


#include <iostream>

NORI_NAMESPACE_BEGIN

class BasicMediumPathIntegrator : public Integrator {
public:
    BasicMediumPathIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        MediumSamplingRecord mRec;
        Color3f tp(1.0f, 1.0f, 1.0f);
        Color3f Li(0.0f, 0.0f, 0.0f);
        Ray3f pathRay(ray.o, ray.d);
        float eta = 1.0f;


        //check for intersection
        bool validIntersection = scene->rayIntersect(ray, its);

        while(true) {
            if(validIntersection && its.mesh->hasMedium()){
                const Medium *curMedium = its.mesh->getMedium();
                Frame w2v(pathRay.d);
                PhaseFunctionQueryRecord pRec;
                bool scatterEvent = curMedium->sampleDistance(Ray3f(pathRay, 0.0, its.t), mRec, sampler->next2D());
                if(scatterEvent) {
                    const PhaseFunction *phase = curMedium->getPhaseFunction();
                    tp *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;


                    pRec = PhaseFunctionQueryRecord(w2v.toLocal(-pathRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);
                    float phaseVal = phase->sample(pRec, sampler->next2D());

                    if (phaseVal == 0){
                        std::cout << "phaseVal is zero!" << std::endl;
                        break;
                    }
                    tp *= phaseVal;

                    //std::cerr << "scattering tp = " << tp << std::endl;
                } else {
                    // no scattering direct fly trough the medium
                    //std::cerr << "Fly trough" << std::endl;
                    tp *= mRec.transmittance / mRec.pdfFailure;
                }
                if(scatterEvent) {
                    Vector3f wo = w2v.toWorld(pRec.wo);
                    pathRay = Ray3f(mRec.p, wo);
                } else {
                    pathRay = Ray3f(mRec.p, pathRay.d);
                }
                validIntersection = scene->rayIntersect(pathRay, its);

            } else if(!validIntersection) {
                const Emitter* distantsDisk = scene->getDistantEmitter();
                //sample the distant disk light
                if(distantsDisk != nullptr ) {
                    tp *= distantsDisk->sampleL(pathRay.d);
                }
                break;
            } else {
                // surface interaction
                //get the radiance of hitten object
                if (its.mesh->isEmitter() ) {
                    const Emitter* areaLightEM = its.mesh->getEmitter();
                    const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
                    Li += tp * aEM->sampleL(-pathRay.d, its.shFrame.n, its);
                }

                //get the asigned BSDF
                const BSDF* curBSDF = its.mesh->getBSDF();

                //transform to the local frame
                BSDFQueryRecord bRec = BSDFQueryRecord(its.toLocal(-pathRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);

                //sample the BSDF
                Color3f bsdfVal =  curBSDF->sample(bRec, sampler->next2D());
                if (bsdfVal.maxCoeff() <= 0.0f)
                    break;

                const Vector3f wo = its.toWorld(bRec.wo);

                float woDotGeoN = its.geoFrame.n.dot(wo);
                if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0.0)
                    break;

                tp *= bsdfVal;

                eta *= bRec.eta;
                pathRay = Ray3f(its.p, wo);
                validIntersection = scene->rayIntersect(pathRay, its);
            }
            // Russian roulette
            float q = std::min(tp.maxCoeff() * eta * eta, (float) 0.95f);
            if (sampler->next1D() >= q)
                break;
            tp /= q;
        }
        return Li;
    }

    std::string toString() const {
        return "BasicMediumPathIntegrator[]";
    }

};

NORI_REGISTER_CLASS(BasicMediumPathIntegrator, "basicMediumPathIntegrator");
NORI_NAMESPACE_END
