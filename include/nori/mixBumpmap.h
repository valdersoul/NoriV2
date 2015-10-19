#include <nori/bumpTexture.h>
#include <nori/perlinTexture.h>
#include <nori/bitmap.h>
#include <filesystem/resolver.h>

NORI_NAMESPACE_BEGIN


class MixBumpmap : public BumpTexture {
public:
    MixBumpmap(const PropertyList &propList);

    ~MixBumpmap() {
        delete m_bumpTexture1;
        delete m_bumpTexture2;
        delete m_perlinNoise;
    }

    Color3f sample(Point2f &uv) const { return Color3f(0.0f);}

    void getNormal(Point2f &uv, Normal3f &n) const;
    void getNormalTangentSpace(Point2f &uv, Normal3f &n) const;

    void addChild(NoriObject *obj);

    std::string toString() const;

    EClassType getClassType() const { return EMIXBUMPMAP; }
private:
    BumpTexture *m_bumpTexture1 = nullptr;
    BumpTexture *m_bumpTexture2 = nullptr;
    PerlinTexture *m_perlinNoise = nullptr;
    float m_threshold = 0.0;

};

NORI_NAMESPACE_END
