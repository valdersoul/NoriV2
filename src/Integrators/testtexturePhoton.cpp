/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>
#include <nori/areaLight.h>

NORI_NAMESPACE_BEGIN

class testTexturePhoton : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    testTexturePhoton(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
        m_shootedRays = 0;
    }

    void preprocess(const Scene *scene) {
        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

        /* Estimate a default photon radius */
        if (m_photonRadius == 0)
            m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

        int storedPhotons = 0;

        const  std::vector<Emitter *> lights = scene->getEmitters();
        int nLights = lights.size();
        Color3f tp(1.0f, 1.0f, 1.0f);

        cout << "Starting to create "<< m_photonCount << " photons!" << endl;
        int percentDone= 0;
        int onePercent = int(floor(m_photonCount / 100.0));

        // create the expected number of photons
        while(storedPhotons < m_photonCount) {
            //uniformly sample 1 light (assuming that we only have area lights)
            int var = int(std::min(sampler->next1D()*nLights, float(nLights) - 1.0f));
            const areaLight* curLight = static_cast<const areaLight *> (lights[var]);

            //sample a photon
            Photon curPhoton;
            Vector3f unQuantDir(0.0f,0.0f,0.0f);
            curLight->samplePhoton(sampler, curPhoton, 1, nLights, unQuantDir);
            Color3f alpha = curPhoton.getPower();
            Color3f tp(1.0f, 1.0f, 1.0f);


            //trace the photon
            Intersection its;
            Ray3f photonRay(curPhoton.getPosition(), unQuantDir);
            m_shootedRays++;
            if (scene->rayIntersect(photonRay, its)) {
                while(true) {
                    const BSDF* curBSDF = its.mesh->getBSDF();


                    if (curBSDF->isDiffuse()) {
                        //store the photon
                        m_photonMap->push_back(Photon(
                            its.p  /* Position */,
                            -photonRay.d /* Direction*/,
                            tp * alpha  /* Power */
                        ));
                        storedPhotons++;
                    }

                    if(!(storedPhotons < m_photonCount)) break;

                    BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-photonRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);
                    Color3f fi =  curBSDF->sample(query, sampler->next2D());
                    Color3f tex(0.0f, 0.0f, 0.0f);
                    if(its.mesh->hasTexture()){
                        if((!isnan(its.uv(0))) && (!isnan(its.uv(1))))
                            tex = its.mesh->getTexture()->sample(its.uv);
                    } else {
                        tex = fi;
                    }

                    if(tex.maxCoeff() == 0.0f) break;

                    tp *= fi;

                    Vector3f wo = its.toWorld(query.wo);
                    photonRay = Ray3f(its.p, wo);

                    //ray escapes the scene
                    if (!scene->rayIntersect(photonRay, its)) break;

                    //stop critirium russian roulette
                    float q = tp.maxCoeff();
                    if(q < sampler->next1D()) break;
                    tp /= q;
                }

            }
            if(onePercent != 0) {
                if(storedPhotons % onePercent == 0){
                    int percent = int(floor(storedPhotons / onePercent));
                    if(percent % 10 == 0 && percentDone != percent){
                        percentDone = percent;
                        cout << percent << "%" << endl;
                    }
                }
            }
        }

        /* Build the photon map */
        m_photonMap->build();
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;


        //check if the ray intersects the scene
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0f);
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

            //Normal3f n = its.shFrame.n;



            if(curBSDF->isDiffuse()) {
                std::vector<uint32_t> results;
                m_photonMap->search(its.p, m_photonRadius, results);
                Color3f Li(0.0f, 0.0f, 0.0f);
                int k = results.size();
                //cout << k << endl;

                if(k > 0) {

                    //cout << results.size() << " Photons found!" << endl;
                    //get the power from all photons
                    //for (uint32_t i : results)
                    //const Photon &photonk = (*m_photonMap)[k-1];
                    Color3f Lindir(0.0f, 0.0f, 0.0f);
                    for (int i = 0; i < k; ++i)
                    {
                        const Photon &photon = (*m_photonMap)[results[i]];
                        Vector3f wi = its.toLocal(photon.getDirection());
                        Vector3f wo = its.toLocal(its.shFrame.n);
                        BSDFQueryRecord dummy = BSDFQueryRecord(wi, wo, EMeasure::ESolidAngle);
                        Color3f f = curBSDF->eval(dummy);

                        Color3f tex(0.0f, 0.0f, 0.0f);
                        if(its.mesh->hasTexture()){
                            if((!isnan(its.uv(0))) && (!isnan(its.uv(1))))
                                tex = its.mesh->getTexture()->sample(its.uv);
                        } else {
                            tex = f;
                        }

                        Lindir += (tp * tex) * photon.getPower() / (M_PI * m_photonRadius * m_photonRadius);


                    }

                    Li += Lindir;

                    //cout << "Li = " <<  Li.toString() << endl;
                    if(Li.maxCoeff() > 0.0f)
                        L +=  Li  / m_shootedRays;

                }
                break;
            }

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

            //ray escapes the scene
            if (!scene->rayIntersect(pathRay, its)) break;


            //stop critirium russian roulette
            float maxCoeff = tp.maxCoeff();
            float q = std::min(0.99f, maxCoeff);
            if(q < sampler->next1D()) break;
            tp /= q;
        }
        return L;
    }

    std::string toString() const {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    int m_photonCount;
    float m_photonRadius;
    int m_shootedRays;
    std::unique_ptr<PhotonMap> m_photonMap;
};

NORI_REGISTER_CLASS(testTexturePhoton, "testTexturePhoton");
NORI_NAMESPACE_END
