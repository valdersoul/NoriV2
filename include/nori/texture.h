#if !defined(__NORI_TEXTURE_H)
#define __NORI_TEXTURE_H

#include <nori/object.h>



NORI_NAMESPACE_BEGIN

class Texture : public NoriObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return ETEXTURE; }

    virtual Color3f sample(Point2f &uv) const = 0;

};

NORI_NAMESPACE_END
#endif /* __NORI_TEXTURE_H */
