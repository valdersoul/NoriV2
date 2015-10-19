#include <nori/emitter.h>
#include <nori/scene.h>
#include <iostream>

NORI_NAMESPACE_BEGIN

class areaLight : public Emitter {
public:
    areaLight (const PropertyList &props);

    virtual std::string toString() const;

    virtual Color3f Power(const Scene *) const;

    virtual bool isDeltaLight() const;


    virtual Color3f sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /*float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const;

    //TEMPORY
    virtual Color3f sampleL(const Vector3f &d, const Normal3f &n, const Intersection &its) const;

    //TEMPORY
    virtual Color3f sampleL(const Vector3f &d) const;

    virtual void setMesh(Mesh* mesh);

    virtual Color3f radiance() const;

    virtual float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const;

    virtual void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const;

protected:
    Color3f m_radiance;
    bool m_shootInNormalDirc;
    Mesh* m_mesh;
};

NORI_NAMESPACE_END
