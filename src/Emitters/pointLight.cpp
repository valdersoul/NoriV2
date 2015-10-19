#include <nori/emitter.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class pointLight : public Emitter {
public:
    pointLight(const PropertyList &props) :
        Emitter()
    {
        //set the arguments
        m_position =  props.getPoint("position");
        m_power = props.getColor("power");
    }

    std::string toString() const {
        return "pointLight[]";
    }

    Color3f Power(const Scene *) const {return  4.f * M_PI * m_power; }

    bool isDeltaLight() const { return true; }


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

        //init the visibility tester
        vis->SetSegment(p, pEpsilon, m_position, 0.0f);

        // calc dist squared between m_position and p
        float distSquared = (m_position - p).squaredNorm();

        //retrun the power divided squared distance and four pi
        return m_power / (distSquared * 4 * M_PI);
    }



    Color3f sampleL(const Vector3f &d) const {
        return Color3f(0.0);
    };

    Color3f radiance() const {
        return m_power * INV_FOURPI;
    };

    float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const {
        return 0.0f;
    }
    void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const {

        //get a direction
        Vector3f dir = Warp::squareToUniformSphere(sampler->next2D());

        unQuantDir = dir;

        dir = LightToWorld * dir;

        ////const Color3f power = (M_PI * A * m_radiance) / (float(N) * float(nLights));
        const Color3f power = (m_power * float(nLights)) / (float(N));

        // create the photon
        photon = Photon(m_position, dir, power);
    }

private:
    Point3f m_position;
    Color3f m_power;
};

NORI_REGISTER_CLASS(pointLight, "point");
NORI_NAMESPACE_END
