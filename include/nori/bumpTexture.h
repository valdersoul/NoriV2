#if !defined(__NORI_BUMPTEXTURE_H)
#define __NORI_BUMPTEXTURE_H

#include <nori/texture.h>



NORI_NAMESPACE_BEGIN


class BumpTexture : public Texture {
public:

    virtual void getNormal(Point2f &uv, Normal3f &n) const = 0;

    virtual void getNormalTangentSpace(Point2f &uv, Normal3f &n) const = 0;

    EClassType getClassType() const { return EBUMPMAP; }


};

NORI_NAMESPACE_END
#endif /* __NORI_BUMPTEXTURE_H */
