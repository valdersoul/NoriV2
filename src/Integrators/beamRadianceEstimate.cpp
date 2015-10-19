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
#include <nori/volume.h>
#include <fstream>

NORI_NAMESPACE_BEGIN

class BeamRadianceEstimate : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    BeamRadianceEstimate(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 100000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
        m_shootedRays = 0;
        m_stepsPerRay = props.getInteger("nrStepMAX", 10 /* Default: automatic */);
        m_k = props.getInteger("k", 10 /* Default: automatic */);
        m_stepSize = props.getFloat("stepSize", 0.1 /* Default: automatic */);
        debug = props.getBoolean("debug", false);
        m_singleScatter = props.getBoolean("singleScatter", false);
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

    void preprocess(const Scene *scene) {
        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMapVolume = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMapSurfaces = std::unique_ptr<PhotonMap>(new PhotonMap());

        m_photonMapVolume->reserve(m_photonCount);
        m_photonMapSurfaces->reserve(m_photonCount);
        const std::vector<Volume *> volumes = scene->getVolumes();

        /* Estimate a default photon radius */
        if (m_photonRadius == 0)
            m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

        int storedPhotons = 0;

        const  std::vector<Emitter *> lights = scene->getEmitters();
        int nLights = lights.size();

        cout << "Starting to create "<< m_photonCount << " photons!" << endl;
        int percentDone= 0;
        int onePercent = int(floor(m_photonCount / 100.0));

        // create the expected number of photons
        while(storedPhotons < m_photonCount) {
            //<CREATE PHOTON>
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

            float pdfMedium;
            float pdfSurface;
            //</CREATE PHOTON>

            //<TRACE PHOTON>
             while(true) {
                bool surIntersection = scene->rayIntersect(photonRay, its);
                //check if there is a medium between the photon and the current surface
                Volume *curVol = nullptr;
                float ffDist = 0.0f;
                float t0 = 0.0f;
                float t1 = 0.0f;
                if(checkVolIntersetcion(photonRay, t0, t1, volumes, curVol)) {
                    int sigIndex = int(std::min(sampler->next1D() * 3.0f, 2.0f));
                    ffDist = curVol->getFreeFlightDistance(sampler->next1D(), sigIndex, pdfMedium, pdfSurface);
                }
                float volEnd = t1 - t0;
                if (curVol != nullptr &&   /* check if there is a volume && t0 <= ffDist */
                    ffDist <= volEnd &&      /* check if the free flight distance is inside the medium */
                    !m_insideMesh &&
                    (!surIntersection || /* no mesh in that direction */
                    ffDist < (its.t - t0))) //  ||  /* mesh inside th volume, but freeflight distance smaller */
                    //t1 < (its.t - t0)))        /* intersection is after the volume */
                {
                    //Medium

                    //compute the postion
                    Point3f o = photonRay(ffDist + t0);
                    //Create a frame to transform to volume
                    Frame w2v(photonRay.d);
                    //the current phase function
                    const PhaseFunction *curPhase = curVol->getPhaseFkt();
                    //cout << tp * alpha << endl;
                    //store the photon
                    m_photonMapVolume->push_back(Photon(
                        o  /* Position */,
                        -photonRay.d /* Direction*/,
                        tp * alpha  /* Power */
                    ));
                    storedPhotons++;
                    if(!(storedPhotons < m_photonCount)) break;

                    PhaseFunctionQueryRecord query(w2v.toLocal(-photonRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);
                    curPhase->sample(query, sampler->next2D());
                    Color3f albido = curVol->getSimga_s() / curVol->getSimga_t();

                    if(albido.maxCoeff() == 0.0f) break;
                    if(m_singleScatter) break;
                    tp *= albido / query.pdf;

                    Vector3f wo = w2v.toWorld(query.wo);
                    photonRay = Ray3f(o, wo);


                } else  {
                    if(!surIntersection) break;

                    const BSDF* curBSDF = its.mesh->getBSDF();

                    if (curBSDF->isDiffuse()) {
                        //store the photon
                        m_photonMapSurfaces->push_back(Photon(
                        its.p  /* Position */,
                        -photonRay.d /* Direction*/,
                        tp * alpha  /* Power */
                        ));
                        storedPhotons++;
                    }

                    if(!(storedPhotons < m_photonCount)) break;

                    BSDFQueryRecord query = BSDFQueryRecord(its.toLocal(-photonRay.d), Vector3f(0.0f), EMeasure::ESolidAngle);
                    Color3f fi =  curBSDF->sample(query, sampler->next2D());

                    if(its.mesh->hasTexture()){
                        if((!isnan(its.uv(0))) && (!isnan(its.uv(1))))
                            fi = its.mesh->getTexture()->sample(its.uv);
                    }

                    m_insideMesh = query.entering;

                    if(fi.maxCoeff() == 0.0f) break;

                    tp *= fi;

                    Vector3f wo = its.toWorld(query.wo);
                    photonRay = Ray3f(its.p, wo);


                }
                //stop critirium russian roulette
                float q = tp.maxCoeff();
                if(q < sampler->next1D()) break;
                    tp /= q;

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
            //</TRACE PHOTON>
        }
        cout << m_photonMapSurfaces->size() << " Photons landed on the surface!" << endl;
        cout << m_photonMapVolume->size() << " Photons landed in the volume!" << endl;


        /* Build the photon maps */
        if(m_photonMapSurfaces->size() > 0)
            m_photonMapSurfaces->build(true);
        if(m_photonMapVolume->size() > 0) {
            m_photonMapVolume->build(true);
            m_photonMapVolume->estimateRadius(size_t(m_k));
            if(debug) {
                std::ofstream of;
                of.open("beamPhotons.off");
                of << "NOFF" << endl;
                of << m_photonMapVolume->size() << " " << 0 << " " << 0 << endl;
                for (unsigned int i = 0; i < m_photonMapVolume->size(); i++){
                    const Photon& ph = (*m_photonMapVolume)[i];
                    Point3f pos = ph.getPosition();
                    Vector3f dir = ph.getDirection();
                    of << tfm::format("%f %f %f %f %f %f", pos.x(), pos.y(), pos.z(),
                        dir.x(), dir.y(), dir.z()) << endl;
                }
                of.close();
            }
            m_photonMapVolume->constructBBH();
            if(debug) {
                BoundingBox3f TEST(Point3f(0.0f, 0.0f, 0.0f), Point3f(0.0f, 0.0f, 0.0f));
                float maxRad = 0.0f;
                float inf = std::numeric_limits<float>::infinity();
                Point3f BBmin(inf, inf, inf);
                Point3f BBmax(-inf, - inf, - inf);

                for (unsigned int i = 0; i < m_photonMapVolume->size(); i++){
                    const Photon& ph = (*m_photonMapVolume)[i];
                    const BoundingBox3f BB = ph.getBB();
                    TEST.expandBy(BB);
                    for (int i = 0; i < 3; ++i) {
                        maxRad = std::max(ph.getRadius(), maxRad);
                        BBmin(i) = std::min(BBmin(i), BB.min(i));
                        BBmax(i) = std::max(BBmax(i), BB.max(i));
                    }
                }
                cout << tfm::format("%f %f %f %f %f %f", BBmin.x(), BBmin.y(), BBmin.z(),
                                        BBmax.x(), BBmax.y(), BBmax.z()) << endl;
                cout << maxRad << endl;
                cout << "--------------"<< endl;
                cout << tfm::format("%f %f %f %f %f %f", TEST.min.x(), TEST.min.y(), TEST.min.z(),
                                        TEST.max.x(), TEST.max.y(), TEST.max.z()) << endl;
            }
        }
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        const std::vector<Volume *> volumes = scene->getVolumes();


        Color3f tp(1.0f, 1.0f, 1.0f);
        Color3f L(0.0f, 0.0f, 0.0f);

        float t0 = 0.0f;
        float t1 = 0.0f;
        Volume *curVol = nullptr;
        Ray3f pathRay(ray.o, ray.d);
        while(true) {
            Color3f Tr(1.0f, 1.0f, 1.0f);
            bool surIntersection = scene->rayIntersect(pathRay, its);
            bool mediumIntersetion = checkVolIntersetcion(pathRay, t0, t1, volumes, curVol);

            //compute the raymarching inside the medium
            if(mediumIntersetion) {
                //Create a frame to transform to volume
                Frame w2v(pathRay.d);

                const PhaseFunction *curPhase = curVol->getPhaseFkt();


                // get the new ray postion
                Point3f mediumBegin = pathRay(t0);
                std::vector<uint32_t> results;
                Ray3f volRay = Ray3f(pathRay.o, pathRay.d, pathRay.mint, t1);
                m_photonMapVolume->gatherPhotons(volRay, results);

                Color3f Li(0.0f, 0.0f, 0.0f);
                int k = results.size();
                std::ofstream of;
                if(debug) {
                    of.open("qryPhotons.off");
                    of << "NOFF" << endl;
                    of << k << " " << 0 << " " << 0 << endl;


                }
                // get the contributions from the photons
                for (int i = 0; i < k; ++i) {
                    const Photon &photon = (*m_photonMapVolume)[results[i]];
                    Vector3f wi = w2v.toLocal(photon.getDirection());
                    Vector3f wo = w2v.toLocal(-pathRay.d);
                    PhaseFunctionQueryRecord pQry(wi, wo, EMeasure::ESolidAngle);
                    curPhase->eval(pQry);

                    float radius = photon.getRadius();
                    float disk = M_PI * radius * radius;

                    Color3f phi = photon.getPower();

                    if(debug)
                        of << tfm::format("%f %f %f %f %f %f", photon.getPosition().x(), photon.getPosition().y(), photon.getPosition().z(),
                            photon.getDirection().x(), photon.getDirection().y(), photon.getDirection().z()) << endl;
                    //get the correct distance
                    Vector3f eyePhoton = photon.getPosition() - mediumBegin;
                    float costheta = std::cos(eyePhoton.dot(volRay.d));
                    float hypo = eyePhoton.norm();
                    float dist = costheta * hypo;

                    Color3f Ti = curVol->getTransmittance(dist);
                    Color3f ss = curVol->sigma_s(photon.getPosition(), - volRay.d);

                    Li += (Ti * ss * pQry.pdf * phi) / disk;
                }

                if(debug){
                    of.close();
                    throw NoriException("BLA");
                }
                Li *= tp;

                if(Li.maxCoeff() > 0.0f)
                    L +=  Li  / m_shootedRays;
                if(!surIntersection)
                    break;
            }
            if(mediumIntersetion) {
                Tr = curVol->getTransmittance(t1 - t0);
                tp *= Tr;
            }

            //get the radiance of hitten object
            if(surIntersection) {
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
                    /*
                    std::vector<uint32_t> results;
                    m_photonMapSurfaces->search(its.p, m_photonRadius, results);
                    Color3f Li(0.0f, 0.0f, 0.0f);
                    int k = results.size();

                    if(k > 0) {
                        Color3f Lindir(0.0f, 0.0f, 0.0f);
                        for (int i = 0; i < k; ++i)
                        {
                            const Photon &photon = (*m_photonMapSurfaces)[results[i]];
                            Vector3f wi = its.toLocal(photon.getDirection());
                            Vector3f wo = its.toLocal(its.shFrame.n);
                            BSDFQueryRecord dummy = BSDFQueryRecord(wi, wo, EMeasure::ESolidAngle);
                            Color3f f = curBSDF->eval(dummy);

                            Lindir += (tp * f) * photon.getPower() / (M_PI * m_photonRadius * m_photonRadius);

                        }

                        Li += Lindir;

                        if(Li.maxCoeff() > 0.0f)
                            L +=  Li  / m_shootedRays;


                    }
                    break;
                    */
                    PhotonMap::SearchResult *results = new PhotonMap::SearchResult[m_k + 1];
                    int k = m_photonMapSurfaces->nnSearch(its.p, m_k, results);

                    Color3f Li(0.0f, 0.0f, 0.0f);

                    if(k > 0) {
                        Color3f Lindir(0.0f, 0.0f, 0.0f);
                        float radius = 0.0f;
                        for (int i = 0; i < k; ++i)
                        {

                            const Photon &photon = (*m_photonMapSurfaces)[results[i].index];
                            Vector3f wi = its.toLocal(photon.getDirection());
                            Vector3f wo = its.toLocal(its.shFrame.n);
                            BSDFQueryRecord dummy = BSDFQueryRecord(wi, wo, EMeasure::ESolidAngle);
                            Color3f f = curBSDF->eval(dummy);
                            if(its.mesh->hasTexture()){
                                if((!isnan(its.uv(0))) && (!isnan(its.uv(1))))
                                    f = its.mesh->getTexture()->sample(its.uv);
                            }

                            radius = std::max((its.p - photon.getPosition()).norm(), radius);
                            Lindir += (tp * f) * photon.getPower();

                        }

                        Li += Lindir / (M_PI * radius * radius);

                        if(Li.maxCoeff() > 0.0f)
                            L +=  Li  / m_shootedRays;

                    }
                    delete results;
                    break;
                }

                //sample the BRDF
                Color3f fi =  curBSDF->sample(query, sampler->next2D());

                if(its.mesh->hasTexture()){
                    if((!isnan(its.uv(0))) && (!isnan(its.uv(1))))
                        fi *= its.mesh->getTexture()->sample(its.uv);
                }

                //check for black brdf
                if(fi.maxCoeff() > 0.0f) {
                    tp *= fi;
                } else {
                    break;
                }
                Vector3f wo = its.toWorld(query.wo);
                pathRay = Ray3f(its.p, wo);

                //stop critirium russian roulette
                float maxCoeff = tp.maxCoeff();
                float q = std::min(0.99f, maxCoeff);
                if(q < sampler->next1D()) break;
                tp /= q;

            } else {
                    break;
            }


        }
        return L;
    }

    std::string toString() const {
        return tfm::format(
            "BeamRadianceEstimate[\n"
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
    int m_stepsPerRay;
    int m_k;
    float m_stepSize;
    std::unique_ptr<PhotonMap> m_photonMapSurfaces;
    std::unique_ptr<PhotonMap> m_photonMapVolume;
    bool debug;
    bool m_insideMesh = false;
    bool m_singleScatter;
};

NORI_REGISTER_CLASS(BeamRadianceEstimate, "beamRadianceEstimate");
NORI_NAMESPACE_END
