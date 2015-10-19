#include <nori/texture.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN


class PerlinTexture : public Texture {
public:
    PerlinTexture(const PropertyList &propList);

    ~PerlinTexture(){
        delete m_sampler;
    }

    Color3f sample(Point2f &uv) const;

    void getNormal(Point2f &uv, Normal3f &n) const;

    std::string toString() const;

    void bilinearInterpolated(float x, float y, float &result) const;

    float turbulence(float x, float y) const;

    void createTexture();

    void addChild(NoriObject *obj);

    EClassType getClassType() const { return EPERLIN; }


private:
    int m_width = 0;
    int m_height = 0;
    int m_levels = 0;
    Sampler *m_sampler = nullptr;
    Eigen::MatrixXf m_texture;

};

NORI_NAMESPACE_END
