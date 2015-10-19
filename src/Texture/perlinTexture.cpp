#include <nori/perlinTexture.h>


NORI_NAMESPACE_BEGIN

    PerlinTexture::PerlinTexture(const PropertyList &propList) {
        m_width = propList.getInteger("width", 256);
        m_height = propList.getInteger("height", 128);
        m_levels = propList.getInteger("levels", 1);

        m_texture.resize(m_height, m_width);


    }
    float PerlinTexture::turbulence(float x, float y) const{
        float value = 0.0;
        float size = m_levels;

        while(size >= 1)
        {
            float curValue = 0.0f;
            bilinearInterpolated(x / size, y / size, curValue);
            value += curValue * size;
            size /= 2.0;
        }

        return (value / m_levels);
    }

    void PerlinTexture::createTexture(){
        for (int x = 0; x < m_width; ++x) {
            for (int y = 0; y < m_height; ++y) {
                m_texture(y, x) = m_sampler->next1D();
            }
        }
    }



    Color3f PerlinTexture::sample(Point2f &uv) const {
        float result = 0.0f;
        float x =  clamp(uv(0), 0.0f, 1.0f) * m_width;
        float y =  clamp(uv(1), 0.0f, 1.0f) * m_height;

        result = turbulence(x, y);
        return Color3f(result, result, result);
    }

    void PerlinTexture::bilinearInterpolated(float x, float y, float &result) const {

        int x1 = int(x);
        int y1 = int(y);
        int x2 = (x1 + 1);
        int y2 = (y1 + 1);
        float Q11 = 0.0f;
        float Q21 = 0.0f;
        float Q12 = 0.0f;
        float Q22 = 0.0f;

        if(x1 < m_width && y1 < m_height){
            Q11 = m_texture(y1, x1);
        }

        if(x2 < m_width && y1 < m_height){
            Q21 = m_texture(y1, x2);;
        }

        if(x1 < m_width && y2 < m_height){
            Q12 = m_texture(y2, x1);
        }

        if(x2 < m_width && y1 < m_height){
            Q22 = m_texture(y2, x2);;
        }

        float fact1 = ((x2 - x) / (x2 - x1));
        float fact2 = ((x - x1) / (x2 - x1));

        float R1 = fact1 * Q11 + fact2 * Q21;
        float R2 = fact1 * Q12 + fact2 * Q22;

        fact1 = ((y2 - y) / (y2 - y1));
        fact2 = ((y - y1) / (y2 - y1));

        result = fact1 * R1 + fact2 * R2;
    }

    void PerlinTexture::addChild(NoriObject *obj){
        switch (obj->getClassType()) {
            case ESampler:
                if (m_sampler)
                    throw NoriException("There can only be one sampler per scene!");
                m_sampler = static_cast<Sampler *>(obj);
                createTexture();
                break;

            default:
                throw NoriException("Scene::addChild(<%s>) is not supported!",
                    classTypeName(obj->getClassType()));
        }
    }


    std::string PerlinTexture::toString() const {
        return "perlinTexture[]";
    }


NORI_REGISTER_CLASS(PerlinTexture, "perlinTexture");
NORI_NAMESPACE_END
