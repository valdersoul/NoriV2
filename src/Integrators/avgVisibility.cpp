#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AverageVisibility : public Integrator {
public:
    AverageVisibility(const PropertyList &props) {
        /* No parameters this time */
        m_length = props.getFloat("length");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(1.0f);

        //get the Normal at the intersection point
        Normal3f n = its.shFrame.n;
        
        //get a random direction and multiply it by the given length
        Vector3f direction = Warp::sampleUniformHemisphere(sampler, n);

        //create a ray
        Ray3f scdRay(its.p, direction);

        //trace the new ray and check for intersection intersection
        Intersection its2;
        if(scene->rayIntersect(scdRay, its2) && its2.t <= m_length){
            // intersection
            return Color3f(0.0f);
        } else {
            return Color3f(1.0f);
        }
            

        
    }

    std::string toString() const {
        return "AverageVisibility[]";
    }
protected:
    float m_length;
};

NORI_REGISTER_CLASS(AverageVisibility, "av");
NORI_NAMESPACE_END