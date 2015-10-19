#include <nori/bumpTexture.h>
#include <Eigen/Geometry>
#include <nori/bitmap.h>
#include <filesystem/resolver.h>

NORI_NAMESPACE_BEGIN

class SingleBumpMap : public BumpTexture{
public:
    SingleBumpMap(const PropertyList &propList) {
        std::string texturePath = propList.getString("texturePath");
        //check if the path exists
        filesystem::path path = getFileResolver()->resolve(texturePath);
        //try to set it
        try {
            if (path.extension() == "exr") {

                Bitmap curText(path.str());

                m_texture = curText;
                m_width = m_texture.cols();
                m_height = m_texture.rows();
            } else {
                cerr << "Fatal error: unknown file \"" << texturePath
                     << "\", expected an extension of type .exr" << endl;
            }
        } catch (const std::exception &e) {
            cerr << "Fatal error: " << e.what() << endl;
        }
    }

    void bilinearInterpolated(float x,float y, Color3f &result) const {

        int x1 = int(x);
        int y1 = int(y);
        int x2 = (x1 + 1);
        int y2 = (y1 + 1);
        Color3f Q11(0.0f);
        Color3f Q21(0.0f);
        Color3f Q12(0.0f);
        Color3f Q22(0.0f);

        if(x1 < m_width && y1 < m_height){
            Q11 = m_texture(y1, x1);
        }

        if(x2 < m_width && y1 < m_height){
            Q21 = m_texture(y1, x2);;
        }

        if(x1 < m_width && y2 < m_height){
            Q12 = m_texture(y2, x1);
        }

        if(x2 < m_width && y2 < m_height){
            Q22 = m_texture(y2, x2);;
        }

        float fact1 = ((x2 - x) / (x2 - x1));
        float fact2 = ((x - x1) / (x2 - x1));

        Color3f R1 = fact1 * Q11 + fact2 * Q21;
        Color3f R2 = fact1 * Q12 + fact2 * Q22;

        fact1 = ((y2 - y) / (y2 - y1));
        fact2 = ((y - y1) / (y2 - y1));

        result = fact1 * R1 + fact2 * R2;
    }

    Color3f sample(Point2f &uv) const {
        Color3f result(0.0f, 0.0f, 0.0f);
        float u = uv(0);
        float v = uv(1);

        u = remainderf(u, 1.0f);
        v = remainderf(v, 1.0f);
        if(u < 0.0f) u = 1.0f + u;
        if(v < 0.0f) v = 1.0f + v;

        //flip the v coordiante
        v = 1.0f - v;

        float x =  clamp(u, 0.0f, 1.0f) * m_width;
        float y =  clamp(v, 0.0f, 1.0f) * m_height;

        bilinearInterpolated(x, y, result);



        return result;
    }


    void getNormal(Point2f &uv, Normal3f &n) const {
        //Sample the map
        Color3f normal = sample(uv);
        Vector3f localNormal(0.0f);
        localNormal(0) = 2.0 * normal(0) - 1.0;
        localNormal(1) = 2.0 * normal(1) - 1.0;
        localNormal(2) = 2.0 * normal(2) - 1.0;
        Vector3f gblN = n;
        //get the tangent frame

        Vector3f gblTan;
        Vector3f c1 = gblN.cross(Vector3f(1.0f, 0.0f, 0.0f));
        Vector3f c2 = gblN.cross(Vector3f(0.0f, 1.0f, 0.0f));

        if( c1.norm() > c2.norm()) {
            gblTan = c1;
        } else {
            gblTan = c2;
        }
        Vector3f gblBiTan = gblN.cross(gblTan);
        //normalize the vectors
        gblBiTan.normalize();
        gblBiTan.normalize();

        //create a basis
        Eigen::Matrix3f tangentBasis;
        tangentBasis.row(0) = gblTan;
        tangentBasis.row(1) = gblBiTan;
        tangentBasis.row(2) = gblN;

        //create the inverse basis
        Eigen::Matrix3f toWorld = tangentBasis.inverse().eval();

        //project the normal to world cooridantes
        n = toWorld * localNormal;
        n.normalize();
    }
    void getNormalTangentSpace(Point2f &uv, Normal3f &n) const {
        Color3f normal = sample(uv);

        n(0) = 2.0 * normal(0) - 1.0;
        n(1) = 2.0 * normal(1) - 1.0;
        n(2) = 2.0 * normal(2) - 1.0;

        n(1) *= -1.0f;
        n.normalized();
    }

    std::string toString() const {
        return "singleBumpMap[]";
    }
private:
    Bitmap m_texture;
    int m_width = 0;
    int m_height = 0;

};

NORI_REGISTER_CLASS(SingleBumpMap, "singleBumpMap");
NORI_NAMESPACE_END
