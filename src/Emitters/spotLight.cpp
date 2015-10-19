#include <nori/emitter.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SpotLight : public Emitter {
public:
    SpotLight(const PropertyList &props) :
        Emitter()
    {
        //set the arguments
        m_position =  props.getPoint("position");
        m_intensity = props.getColor("Intensity");
        bool useLookAt = props.getBoolean("useLookAt", true);
        if(useLookAt){
            Vector3f lookPT = props.getPoint("lookAt", Point3f(0.0f, 0.0f, 0.0f));
            m_direction = (lookPT - m_position).normalized();
        } else {
            m_direction = props.getVector("direction");
        }
        m_theta = props.getFloat("theta");
        m_cosFalloffStart = props.getFloat("falloff");

    }

    std::string toString() const {
        return "SpotLight[]";
    }

    Color3f Power(const Scene *) const {return  m_intensity * 2.f * M_PI *  (1.f - .5f * (m_cosFalloffStart + m_theta)); }
    Color3f Power() const {return  m_intensity * 2.f * M_PI *  (1.f - .5f * (m_cosFalloffStart + m_theta)); }

    bool isDeltaLight() const { return true; }

    float Falloff(const Vector3f &w) const {
        Vector3f wl = (WorldToLight * w).normalized();
        float costheta = wl.dot(m_direction);
        if (costheta < m_theta)     return 0.;
        if (costheta > m_cosFalloffStart)   return 1.;
        // Compute falloff inside spotlight cone
        float delta = (costheta - m_theta) /
                      (m_cosFalloffStart - m_theta);
        return delta*delta*delta*delta;
    }


    Color3f sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /*float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const {
        //set the direction to the light
        *wi = (m_position - p).normalized();

        //pdf i one, because it is a point light
        *pdf = 1.f;
        if(wi->dot(m_direction) > m_theta)
            *pdf = 0.0f;

        //init the visibility tester
        vis->SetSegment(p, pEpsilon, m_position, 0.0f);

        // calc dist squared between m_position and p
        float distSquared = (m_position - p).squaredNorm();

        //retrun the power divided squared distance and four pi
        return Falloff(- *wi) * m_intensity / (distSquared);
    }



    Color3f sampleL(const Vector3f &d) const {
        return Color3f(0.0);
    }

    Color3f radiance() const {
        return m_intensity * INV_FOURPI;
    }

    float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const {
        return 0.0f;
    }
    void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const {

        const Point2f ls2 = sampler->next2D();
        const Point3f pos = m_position;

        //get a direction
        const Vector3f cosDir = Warp::squareToUniformSphereCap(ls2, m_theta);
        // transform it to world coordinates
        Frame globalFrame = Frame(m_direction);
        const Vector3f dir = globalFrame.toWorld(cosDir);
        unQuantDir = dir;

        const Color3f power = (Falloff(dir) * Power() * float(nLights)) / (float(N));

        // create the photon
        photon = Photon(pos, dir, power);
    }

private:
    Point3f m_position;
    Vector3f m_direction;
    float m_theta;
    float m_cosFalloffStart;
    Color3f m_intensity;
};

NORI_REGISTER_CLASS(SpotLight, "spotLight");
NORI_NAMESPACE_END
