#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>
#include <nori/volume.h>


#include <iostream>

NORI_NAMESPACE_BEGIN

class VPTWithShadowRays : public Integrator {
public:
    VPTWithShadowRays(const PropertyList &props) {
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
        Color3f tp(1.0f, 1.0f, 1.0f);
        Color3f L(0.0f, 0.0f, 0.0f);
        Ray3f pathRay(ray.o, ray.d);
        const std::vector<Volume *> volumes = scene->getVolumes();
        const std::vector<Emitter *> lights = scene->getEmitters();
        uint32_t nLights = lights.size();
        bool deltaFlag = true;

        while(true) {
            //check if there is a volume ray intersection
            Volume *curVol = nullptr;
            float t0 = 0.0f;
            float t1 = 0.0f;
            bool surIntersection = scene->rayIntersect(pathRay, its);


            //get the free flight distance
            float pdfMedium;
            float pdfSurface;
            float ffDist = 0.0f;
            if(checkVolIntersetcion(pathRay, t0, t1, volumes, curVol)) {
                int sigIndex = int(std::min(sampler->next1D() * 3.0f, 2.0f));
                ffDist = curVol->getFreeFlightDistance(sampler->next1D(), sigIndex, pdfMedium, pdfSurface);
            }
            //Light sampling
            //randomly select a lightsource
            uint32_t var = uint32_t(std::min(sampler->next1D()*nLights, float(nLights) - 1.0f));

            //init the light color
            Color3f Li(0.0f, 0.0f, 0.0f);
            Color3f Ld(1.0f, 1.0f, 1.0f);

            //create a sample for the light
            const Point2f lightSample = sampler->next2D();
            VisibilityTester vis;
            Vector3f wo;
            float lightpdf;
            float bsdfpdf;


            //check if a volume exists, if the free fliegt distance is
            //inside the media and befor the next surface
            if (curVol != nullptr &&   /* check if there is a volume && t0 <= ffDist */
                ffDist <= t1 &&      /* check if the free flight distance is inside the medium */
                (!surIntersection || /* no mesh in that direction */
                 ffDist < (its.t - t0)  ||  /* mesh inside th volume, but freeflight distance smaller */
                 t1 < (its.t - t0)))        /* intersection is after the volume */
            {

                //compute the postion
                Point3f o = pathRay(ffDist + t0);

                //Create a frame to transform to volume
                Frame w2v(pathRay.d);
                //the current phase function
                const PhaseFunction *curPhase = curVol->getPhaseFkt();

                Color3f albido = curVol->getSimga_s();// / curVol->getSimga_t();
                Color3f transmitance = curVol->getTransmittance(t1 - t0);

                if(pdfMedium != 0.0f)
                    tp *= albido *  transmitance  / pdfMedium;

                // <LIGHT SAMPLING>
                Li = lights[var]->sampleL(o, Epsilon, lightSample , &wo, &lightpdf, &vis);               
                lightpdf /= float(nLights);

                // get the probability for this direction
                PhaseFunctionQueryRecord phaseLightQry(w2v.toLocal(-pathRay.d), w2v.toLocal(wo), EMeasure::ESolidAngle);
                float phaseVal = curPhase->eval(phaseLightQry);


                if(lightpdf != 0.0f && phaseVal != 0.0f && vis.Unoccluded(scene)) {                    
                    //compute the transmitance
                    Color3f T(1.0f, 1.0f, 1.0f); //Transmitance
                    Volume *curVolPhase = nullptr;
                    Ray3f ligthRay(o, w2v.toWorld(phaseLightQry.wo), 0.0, std::numeric_limits<float>::infinity());
                    if(checkVolIntersetcion(ligthRay, t0, t1, volumes, curVolPhase)) {
                        T = curVolPhase->getTransmittance(t1 - t0);
                    }

                    float weight = BalanceHeuristic(float(1), phaseLightQry.pdf, float(1), lightpdf);
                    Ld = (T * Li) / lightpdf;
                    L += tp * Ld * weight * phaseLightQry.pdf;
                }
                // </LIGHT SAMPLING>

                // <PHASE SAMPLING>
                //create a query for the phase function
                PhaseFunctionQueryRecord pQry(w2v.toLocal(-pathRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);
                phaseVal = curPhase->sample(pQry, sampler->next2D());
                if(phaseVal != 0.0f){

                    Ray3f shadowRay(o, w2v.toWorld(pQry.wo), 0.0, std::numeric_limits<float>::infinity());
                    Intersection lightIsect;
                    Color3f T(1.0f, 1.0f, 1.0f);

                    // check for media
                    Volume *curVolPhase = nullptr;
                    if(checkVolIntersetcion(shadowRay, t0, t1, volumes, curVolPhase)) {
                        T = curVolPhase->getTransmittance(t1- t0);
                    }

                    lightpdf = 0.0f;
                    if (scene->rayIntersect(shadowRay, lightIsect)) {
                        if(lightIsect.mesh->isEmitter()){
                           const Emitter* areaLightEMcur = lightIsect.mesh->getEmitter();
                           const areaLight* aEMcur = static_cast<const areaLight *> (areaLightEMcur);                           

                           Li = T * aEMcur->sampleL(-shadowRay.d, lightIsect.shFrame.n, lightIsect);
                           lightpdf = aEMcur->pdf(o, (lightIsect.p - o).normalized(), lightIsect.p, Normal3f(lightIsect.shFrame.n));
                        }
                    } else {
                        const Emitter* distantsDisk = scene->getDistantEmitter();
                        if(distantsDisk != nullptr ) {
                            //check if THIS is right!
                            Li = T * distantsDisk->sampleL(lightIsect.toWorld(pQry.wo));
                            lightpdf = distantsDisk->pdf(Point3f(0.0f), wo, Point3f(0.0f), Normal3f(0.0f));
                        }
                    }
                    lightpdf /= float(nLights);

                    float weight = BalanceHeuristic(float(1), pQry.pdf, float(1), lightpdf);

                    if(lightpdf != 0.0f && phaseVal != 0.0f) {
                        //tp *= phaseVal;
                        Ld = weight * Li;
                        L += tp * Ld;

                    }
                } else {
                    break;
                }
                // </PHASE SAMPLING>

                pathRay = Ray3f(o, w2v.toWorld(pQry.wo), 0.0, std::numeric_limits<float>::infinity());

            } else {

                if (curVol != nullptr && pdfSurface != 0.0f){
                    tp *= curVol->getTransmittance(t1 - t0) / pdfSurface;
                }

                if (!surIntersection) {
                    const Emitter* distantsDisk = scene->getDistantEmitter();
                    if(distantsDisk != nullptr ) {
                        //sample the distant disk light
                        Vector3f d = pathRay.d;
                        Color3f T(1.0f, 1.0f, 1.0f);
                        Volume *curVolPhase = nullptr;
                        if(checkVolIntersetcion(pathRay, t0, t1, volumes, curVolPhase)) {
                            T = curVolPhase->getTransmittance(t1- t0);
                        }
                        L += tp * T * distantsDisk->sampleL(d);
                    }
                    break;
                }

                if (its.mesh->isEmitter() && deltaFlag) {
                    const Emitter* areaLightEM = its.mesh->getEmitter();
                    const areaLight* aEM = static_cast<const areaLight *> (areaLightEM);
                    L += tp * aEM->sampleL(-pathRay.d, its.shFrame.n, its);
                }

                const BSDF* curBSDF = its.mesh->getBSDF();
                deltaFlag = curBSDF->isDeltaBSDF();
                //sample the light
                {
                    Li = lights[var]->sampleL(its.p, Epsilon, lightSample , &wo, &lightpdf, &vis);
                    lightpdf /= float(nLights);
                    //check if the pdf of the sample is greater than 0 and if the color is not black
                    if(lightpdf > 0 && Li.maxCoeff() != 0.0f) {

                        //calculate the cosine term wi in my case the vector to the light
                        float cosTerm = std::abs((its.shFrame.n).dot(wo));
                        const BSDFQueryRecord queryEM = BSDFQueryRecord(its.toLocal(- pathRay.d), its.toLocal(wo), EMeasure::ESolidAngle);
                        Color3f f = curBSDF->eval(queryEM);

                        if(f.maxCoeff() != 0.0f && vis.Unoccluded(scene)) {
                            bsdfpdf = curBSDF->pdf(queryEM);
                            float weight = BalanceHeuristic(float(1), lightpdf, float(1), bsdfpdf);
                            if(curBSDF->isDeltaBSDF())  weight = 1.0f;
                            // check for media
                            Color3f T(1.0f, 1.0f, 1.0f); //Transmitance
                            if(checkVolIntersetcion(vis.r, t0, t1, volumes, curVol)) {                                
                                T = curVol->getTransmittance(t1- t0);
                            }

                            if(bsdfpdf > 0.0f) {
                                Ld = (weight * f * T * Li * cosTerm) / lightpdf;
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
                    tp *= fi;
                    if(bsdfpdf > 0.0f) {
                        Ray3f shadowRay(its.p, its.toWorld(queryMats.wo));
                        Intersection lightIsect;

                         if (scene->rayIntersect(shadowRay, lightIsect)) {
                             if(lightIsect.mesh->isEmitter()){
                                const Emitter* areaLightEMcur = lightIsect.mesh->getEmitter();
                                const areaLight* aEMcur = static_cast<const areaLight *> (areaLightEMcur);
                                // check for media
                                Color3f T(1.0f, 1.0f, 1.0f); //Transmitance
                                if(checkVolIntersetcion(shadowRay, t0, t1, volumes, curVol)) {                                   
                                    T = curVol->getTransmittance(t1- t0);
                                }
                                Li = T * aEMcur->sampleL(-shadowRay.d, lightIsect.shFrame.n, lightIsect);
                                lightpdf = aEMcur->pdf(its.p, (lightIsect.p - its.p).normalized(), lightIsect.p, Normal3f(lightIsect.shFrame.n));
                             }
                         } else {
                             const Emitter* distantsDisk = scene->getDistantEmitter();
                             if(distantsDisk != nullptr ) {
                                 //check if THIS is right!
                                 Color3f T(1.0f, 1.0f, 1.0f); //Transmitance
                                 if(checkVolIntersetcion(shadowRay, t0, t1, volumes, curVol)) {
                                     T = curVol->getTransmittance(t1- t0);
                                 }
                                 Li = T * distantsDisk->sampleL(lightIsect.toWorld(queryMats.wo));
                                 lightpdf = distantsDisk->pdf(Point3f(0.0f), wo, Point3f(0.0f), Normal3f(0.0f));
                             }
                         }
                         lightpdf /= float(nLights);
                         //calculate the weights
                         float weight = BalanceHeuristic(float(1), bsdfpdf, float(1), lightpdf);


                         //check if the lightcolor is not black
                         if(Li.maxCoeff() > 0.0f  && lightpdf > 0.0f ) {
                             //Li /= lightpdf; //????
                             //wo in my case the vector to the light
                             Ld = weight * Li;
                             L += tp * Ld;
                         }
                    }



                } else {
                    break;
                }

                wo = its.toWorld(queryMats.wo);
                pathRay = Ray3f(its.p, wo);
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
        return "VPTWithShadowRays[]";
    }

};

NORI_REGISTER_CLASS(VPTWithShadowRays, "vPTWithShadowRays");
NORI_NAMESPACE_END
