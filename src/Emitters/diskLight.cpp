#include <nori/emitter.h>
#include <nori/scene.h>
#include <iostream>

NORI_NAMESPACE_BEGIN

class diskLight : public Emitter {
public:
    diskLight (const PropertyList &props) :
        Emitter(props.getTransform("toWorld", Transform()))
    {
        //set the arguments
        m_radiance = props.getColor("radiance");
        m_theta = props.getFloat("thetaA");

        //convert 2 radians
        m_theta *= M_PI / 180.0;
    }
    //IMPORTANT FUNCTION
    std::string toString() const {
        return "diskLight";
    }

    Color3f Power(const Scene *) const {return  Color3f(0.0f); }

    bool isDeltaLight() const { return m_theta == 0.0? true : false; }


    Color3f sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /*float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const {

        //set the direction to the light
        Vector3f sampleDirection =  Warp::squareToUniformSphereCap(ls,std::cos(m_theta));

        //transform the sample to world coordinates
        *wi = LightToWorld * sampleDirection;


        //set the pdf of the sample
        *pdf = Warp::squareToUniformSphereCapPdf(sampleDirection, std::cos(m_theta));

        //init the visibility tester
        vis->SetRay(p, *wi);
        return 0.0 < *pdf? m_radiance : Color3f(0.0);

    }


    //TEMPORY
    Color3f sampleL(const Vector3f &d) const {

        Vector3f localVector = WorldToLight * d;
        float cosValue = localVector.dot(Vector3f(0.0, 0.0, 1.0));

        if(cosValue >= std::cos(m_theta)) {
            return m_radiance;
        } else {
            return Color3f(0.0f);
        }
    }

    Color3f radiance() const {
        return m_radiance;
    }

    float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const {
        return Warp::squareToUniformSphereCapPdf(WorldToLight * w, std::cos(m_theta));;
    }
    void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const {
        throw NoriException("Not implemented!");
    }



private:
    Color3f m_radiance;
    float m_theta;
};

NORI_REGISTER_CLASS(diskLight , "distantdisk");
NORI_NAMESPACE_END
