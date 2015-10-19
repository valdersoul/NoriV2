#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/areaLight.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class testtexturePath : public Integrator {
public:
    testtexturePath(const PropertyList &props) {

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

                Color3f tp(1.0f, 1.0f, 1.0f);
                Color3f L(0.0f, 0.0f, 0.0f);
                Ray3f pathRay(ray.o, ray.d);
                while(true) {

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
                        //check if the mesh has texture
                        if(its.mesh->hasTexture()){
                            Color3f tex(0.0f);
                            if((!isnan(its.uv(0))) && (!isnan(its.uv(1)))){
                                //get the texture value
                                tex = its.mesh->getTexture()->sample(its.uv);
                                float cosTerm = std::abs(Frame::cosPhi(query.wi));
                                float pdf = curBSDF->pdf(query);
                                if(pdf != 0.0f) fi = (tex * cosTerm) / pdf;
                            }
                        }
                        tp *= fi;
                    } else {
                        //stop
                        // hit a black brdf
                        break;
                    }
                    Vector3f wo = its.toWorld(query.wo);
                    pathRay = Ray3f(its.p, wo);

                    //ray escapes the scene
                    if (!scene->rayIntersect(pathRay, its)) {
                        //cout << "escape" << endl;
                        //check for distant disk light
                        const Emitter* distantsDisk = scene->getDistantEmitter();

                        if(distantsDisk != nullptr ) L += tp * distantsDisk->sampleL(wo);
                        break;
                    }

                    //stop critirium russian roulette
                    float maxCoeff = tp.maxCoeff();
                    float q = std::min(0.99f, maxCoeff);
                    if(q < sampler->next1D()){
                        break;
                    }
                    tp /= q;


                }
                return L;
    }

    std::string toString() const {
        return "testtexturePath []";
    }
private:

};

NORI_REGISTER_CLASS(testtexturePath, "testtexturePath");
NORI_NAMESPACE_END
