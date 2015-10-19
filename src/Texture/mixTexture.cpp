#include <nori/mixTexture.h>

NORI_NAMESPACE_BEGIN



    MixTexture::MixTexture(const PropertyList &propList) {
        m_threshold = propList.getFloat("threshold", 0.5);
    }

    Color3f MixTexture::sample(Point2f &uv) const {

        Color3f p = m_perlinNoise->sample(uv);
        float perlinSample = p(0);


        if(perlinSample < m_threshold){
            return m_texture1->sample(uv);
        } else {
            return m_texture2->sample(uv);
        }
    }


    void MixTexture::addChild(NoriObject *obj){
        switch (obj->getClassType()) {
            case ETEXTURE:
            {
                if (!m_texture1) {
                    m_texture1 = static_cast<Texture *>(obj);
                } else if (!m_texture2) {
                    m_texture2 = static_cast<Texture *>(obj);
                } else {
                    throw NoriException(
                        "Mixtexture: tried to register more than 2 texture instances!");
                }
            }
            break;

            case EPERLIN: {
                    if (m_perlinNoise)
                        throw NoriException(
                            "Mixtexture: tried to register multiple perlin noise instances!");
                    m_perlinNoise = static_cast<PerlinTexture *>(obj);;
                }
                break;

            default:
                throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }


    std::string MixTexture::toString() const {
        return "mixTexture[]";
    }



NORI_REGISTER_CLASS(MixTexture, "mixTexture");
NORI_NAMESPACE_END
