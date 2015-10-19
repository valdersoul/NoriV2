#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/areaLight.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class TextureEmitter : public areaLight {
public:
    TextureEmitter(const PropertyList &props) :
        areaLight(props)
    {
        //set the arguments
        std::string filePath = props.getString("lightProbPath", "");

        filesystem::path path = getFileResolver()->resolve(filePath);
        //try to set it
        try {
            if (path.extension() == "exr") {
                m_lightprobe = new Bitmap(path.str());
                m_wrapper = new HSWrapper(*m_lightprobe);


                m_BBMin = props.getPoint("BBMin");

                m_XVec = props.getVector("XVec", Vector3f(1.0f, 0.0f, 0.0f));
                m_YVec = props.getVector("YVec", Vector3f(0.0f, 1.0f, 0.0f));
                bool flipNormal = props.getBoolean("flipNormal", true);
                m_intensity = props.getColor("Intensity", Color3f(1.0f, 1.0f, 1.0f));
                if(props.getBoolean("computeNormal", true)) {
                    m_ns = m_XVec.cross(m_YVec);
                } else {
                    m_ns = props.getVector("normal", Vector3f(0.0f, 0.0f, -1.0f));
                }
                if(flipNormal)
                    m_ns = -m_ns;
                m_ns.normalize();
                m_width = m_lightprobe->cols();
                m_height = m_lightprobe->rows();
            } else {
                cerr << "Fatal error: unknown file \"" << filePath
                     << "\", expected an extension of type .exr" << endl;
            }
        } catch (const std::exception &e) {
            cerr << "Fatal error: " << e.what() << endl;
        }
    }

    std::string toString() const {
        return "TextureEmitter[]";
    }

    Color3f Power(const Scene *) const {return Color3f(0.0f); }

    bool isDeltaLight() const { return false; }


    Color3f sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /*float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const {

        // importance sample a pixel coordiante
        Point2f samplePix =  Warp::squareToHSW(ls, m_wrapper);
        samplePix(1) = 1.0f - samplePix(1);
        Point2f normalizedPix  = samplePix;

        samplePix(0) *= m_width;
        samplePix(1) *= m_height;


        // get the point in world space
        Point3f lightPT = m_BBMin + normalizedPix(0) * m_XVec;
        lightPT += normalizedPix(1) * m_YVec;

        //set the direction to the light
        *wi = (lightPT - p).normalized();
        //init the visibility tester
        vis->SetSegment(p, pEpsilon, lightPT, pEpsilon);

        float imgPDF = Warp::squareToHSWPDF(normalizedPix, m_wrapper);
        //imgPDF *= m_XVec.norm() * m_YVec.norm();
        imgPDF *= m_height * m_width;

        //Check if you need to add the area of the image
        *pdf = imgPDF;

        int x = int(std::min(m_width - 1.0f, samplePix(0)));
        int y = int(std::min(m_height - 1.0f, samplePix(1)));


        Color3f result = m_wrapper->m_lightprob(y, x);
        return 0.0f < *pdf? result : Color3f(0.0);

    }
    ~TextureEmitter(){
        delete m_wrapper;
        delete m_lightprobe;
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
            Q11 = m_wrapper->m_lightprob(y1, x1);
        }

        if(x2 < m_width && y1 < m_height){
            Q21 = m_wrapper->m_lightprob(y1, x2);;
        }

        if(x1 < m_width && y2 < m_height){
            Q12 = m_wrapper->m_lightprob(y2, x1);
        }

        if(x2 < m_width && y2 < m_height){
            Q22 = m_wrapper->m_lightprob(y2, x2);;
        }

        float fact1 = ((x2 - x) / (x2 - x1));
        float fact2 = ((x - x1) / (x2 - x1));

        Color3f R1 = fact1 * Q11 + fact2 * Q21;
        Color3f R2 = fact1 * Q12 + fact2 * Q22;

        fact1 = ((y2 - y) / (y2 - y1));
        fact2 = ((y - y1) / (y2 - y1));

        result = fact1 * R1 + fact2 * R2;
    }


    Color3f sampleL(const Vector3f &d, const Normal3f &n, const Intersection &its) const {
        Color3f result(0.0f, 0.0f, 0.0f);
        if(d.dot(n) > 0 && (!isnan(its.uv(0))) && (!isnan(its.uv(1)))){
            float u = its.uv(0);
            float v = its.uv(1);

            u = remainderf(u, 1.0f);
            v = remainderf(v, 1.0f);
            if(u < 0.0f) u = 1.0f + u;
            if(v < 0.0f) v = 1.0f + v;

            //flip the v coordiante
            v = 1.0f - v;

            float x =  clamp(u, 0.0f, 1.0f) * m_width;
            float y =  clamp(v, 0.0f, 1.0f) * m_height;


            bilinearInterpolated(x, y, result);

        }
        return result;
    }


    Color3f sampleL(const Vector3f &d) const {
        return Color3f(0.0f);
    }
    void setMesh(Mesh* mesh){
        m_mesh = mesh;
    }

    Color3f radiance() const {
        return Color3f(0.0f);
    }

    float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const {
        return 0.0f;
    }
    void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const {

        Point3f ps;
        const Point2f ls1 = sampler->next2D();
        // importance sample a pixel coordiante
        Point2f samplePix =  Warp::squareToHSW(ls1, m_wrapper);
        samplePix(1) = 1.0f - samplePix(1);
        Point2f normalizedPix  = samplePix;

        samplePix(0) *= m_width;
        samplePix(1) *= m_height;


        // get the point in world space
        Point3f lightPT = m_BBMin + normalizedPix(0) * m_XVec;
        lightPT += normalizedPix(1) * m_YVec;

        //get the surface area
        float A = m_XVec.norm() * m_YVec.norm();

        const Point2f ls2 = sampler->next2D();
        const Point3f pos = lightPT;

        //get a direction
        Vector3f cosDir(0.0f, 0.0f, 1.0f);
        if(!m_shootInNormalDirc){
            cosDir = Warp::squareToCosineHemisphere(ls2);
        } else {
            A *= M_PI;
        }
        // transform it to world coordinates
        Frame globalFrame = Frame(m_ns);
        const Vector3f dir = globalFrame.toWorld(cosDir);
        unQuantDir = dir;
        int x = int(std::min(m_width - 1.0f, samplePix(0)));
        int y = int(std::min(m_height - 1.0f, samplePix(1)));

        Color3f radiance = m_wrapper->m_lightprob(y, x);

        const Color3f power = (A * m_intensity * radiance * float(nLights)) / (float(N));

        // create the photon
        photon = Photon(pos, dir, power);
    }


private:
    HSWrapper *m_wrapper = nullptr;
    Bitmap *m_lightprobe = nullptr;
    Point3f m_BBMin;
    Vector3f m_XVec;
    Vector3f m_YVec;
    Vector3f m_ns;
    Color3f m_intensity;
    int m_width = 0;
    int m_height = 0;
};

NORI_REGISTER_CLASS(TextureEmitter, "textureEmitter");
NORI_NAMESPACE_END
