#include <nori/texture.h>
#include <nori/bumpTexture.h>
#include <nori/perlinTexture.h>

NORI_NAMESPACE_BEGIN


class MixTexture : public Texture {
public:
    MixTexture(const PropertyList &propList);

    ~MixTexture() {
        delete m_texture1;
        delete m_texture2;
        delete m_perlinNoise;
    }

    Color3f sample(Point2f &uv) const;


    void addChild(NoriObject *obj);


    std::string toString() const;

    EClassType getClassType() const { return EMIXTEXTURE; }
private:
    Texture *m_texture1 = nullptr;
    Texture *m_texture2 = nullptr;
    PerlinTexture *m_perlinNoise = nullptr;
    float m_threshold = 0.0;


};

NORI_NAMESPACE_END
