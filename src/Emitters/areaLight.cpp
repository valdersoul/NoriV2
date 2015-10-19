#include <nori/areaLight.h>

NORI_NAMESPACE_BEGIN

    areaLight::areaLight (const PropertyList &props) :
        Emitter()
    {
        //set the arguments
        m_radiance = props.getColor("radiance", Color3f(10.0f, 10.0f, 10.0f));
        m_shootInNormalDirc = props.getBoolean("shootInNormal", false);
    }

    std::string areaLight::toString() const {
        return "areaLight";
    }

    Color3f areaLight::Power(const Scene *) const {return  Color3f(0.0f); }

    bool areaLight::isDeltaLight() const { return false; }


    Color3f areaLight::sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /*float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const {

        Normal3f ns;
        Point3f ps;
        //sample the mesh
        m_mesh->samplePosition(ls,ps,ns);

        //set the direction to the light
        *wi = (ps - p).normalized();

        //set the pdf of the sample
        *pdf = m_mesh->Pdf(p, ps, ns, *wi);

        if(isnan(*pdf)){
            cout << "p = " << p << endl;
            cout << "ps = " << ps << endl;
            cout << "ns = " << ns << endl;
            cout << "wi = " << *wi << endl;
        }


        //init the visibility tester
        vis->SetSegment(p, pEpsilon, ps , Epsilon);

        //check if a backfac is hitn
        if(ns.dot(-*wi) < 0.0) return Color3f(0.0, 0.0, 0.0);

        return 0.0 < *pdf? m_radiance : Color3f(0.0, 0.0, 0.0);

    }

    //TEMPORY
    Color3f areaLight::sampleL(const Vector3f &d, const Normal3f &n, const Intersection &its) const {
        return  (d.dot(n) > 0) ? m_radiance : Color3f(0.0f, 0.0f, 0.0f);
    }

    //TEMPORY
    Color3f areaLight::sampleL(const Vector3f &d) const {
            //std::cout << "area hit" << std::endl;
            return m_radiance;
    }

    void areaLight::setMesh(Mesh* mesh){
        m_mesh = mesh;
    }

    Color3f areaLight::radiance() const {
        return m_radiance;
    }
    float areaLight::pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const {
        return m_mesh->Pdf(p, hitPoint, n, w);
    }

    void areaLight::samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const {
        Normal3f ns;
        Point3f ps;
        const Point2f ls1 = sampler->next2D();
        //sample the mesh
        m_mesh->samplePosition(ls1, ps, ns);

        //get the surface area
        float A = m_mesh->m_surfaceArea;

        const Point2f ls2 = sampler->next2D();
        const Point3f pos = ps;

        //get a direction
        Vector3f cosDir(0.0f, 0.0f, 1.0f);
        if(!m_shootInNormalDirc) {
            cosDir = Warp::squareToCosineHemisphere(ls2);
        } else {
            A *= M_PI;
        }
        // transform it to world coordinates
        Frame globalFrame = Frame(ns);
        const Vector3f dir = globalFrame.toWorld(cosDir);
        unQuantDir = dir;

        ////const Color3f power = (M_PI * A * m_radiance) / (float(N) * float(nLights));
        const Color3f power = (A * m_radiance * float(nLights)) / (float(N));

        // create the photon
        photon = Photon(pos, dir, power);
    }





NORI_REGISTER_CLASS(areaLight , "area");
NORI_NAMESPACE_END
