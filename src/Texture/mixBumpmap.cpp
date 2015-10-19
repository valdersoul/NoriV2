#include <nori/mixBumpmap.h>

NORI_NAMESPACE_BEGIN



    MixBumpmap::MixBumpmap(const PropertyList &propList) {
        m_threshold = propList.getFloat("threshold", 0.5);
    }

    void MixBumpmap::getNormal(Point2f &uv, Normal3f &n) const {
        Color3f p = m_perlinNoise->sample(uv);

        if(p(0) < m_threshold){
            return m_bumpTexture1->getNormal(uv, n);
        } else {
            return m_bumpTexture2->getNormal(uv, n);;
        }
    }
    void MixBumpmap::getNormalTangentSpace(Point2f &uv, Normal3f &n) const {
        Color3f p = m_perlinNoise->sample(uv);

        if(p(0) < m_threshold){
            return m_bumpTexture1->getNormalTangentSpace(uv, n);
        } else {
            return m_bumpTexture2->getNormalTangentSpace(uv, n);;
        }
    }

    void MixBumpmap::addChild(NoriObject *obj){
        switch (obj->getClassType()) {
            case EBUMPMAP:{
                if (!m_bumpTexture1) {
                    m_bumpTexture1 = static_cast<BumpTexture *>(obj);
                } else if (!m_bumpTexture2) {
                    m_bumpTexture2 = static_cast<BumpTexture *>(obj);
                } else {
                    throw NoriException(
                        "Mesh: tried to register more than 2 texture instances!");
                }
            }
            break;
            case EPERLIN: {
                    if (m_perlinNoise)
                        throw NoriException(
                            "Mesh: tried to register multiple perlin noise instances!");
                    m_perlinNoise = static_cast<PerlinTexture *>(obj);;
            }
            break;

            default:
                throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }


    std::string MixBumpmap::toString() const {
        return "mixBumpmap[]";
    }



NORI_REGISTER_CLASS(MixBumpmap, "mixBumpmap");
NORI_NAMESPACE_END
