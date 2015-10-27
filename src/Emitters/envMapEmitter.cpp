#include <nori/emitter.h>
#include <nori/scene.h>


NORI_NAMESPACE_BEGIN
using namespace std;

struct LatLongMap {

    LatLongMap(){

    }

    float sign(float x){
        if(x < 0){
            return -1.0f;
        } else {
            return 1.0f;
        }
    }

    Point2f latLong (const Vector3f &dir){
        float r = sqrt (dir(2) * dir(2) + dir(0) * dir(0));

        float latitude = (r < abs (dir(1)))?
                 acos (r / dir.norm()) * sign(dir(1)):
                 asin (dir(1) / dir.norm());

        float longitude = (dir(2) == 0 && dir(0) == 0)? 0: atan2 (dir(0), dir(2));


        return Point2f(latitude, longitude);
    }


    Point2f latLong (const BoundingBox2i &dataWindow, const Vector2f &pixelPosition){
        float latitude, longitude;

        if (dataWindow.max(1) > dataWindow.min(1)){
        latitude = -M_PI *
              ((pixelPosition(1)  - dataWindow.min(1)) /
               (dataWindow.max(1) - dataWindow.min(1)) - 0.5f);
        } else {
            latitude = 0;
        }

        if (dataWindow.max(0) > dataWindow.min(0)) {
        longitude = -2 * M_PI *
               ((pixelPosition(0)  - dataWindow.min(0)) /
                (dataWindow.max(0) - dataWindow.min(0)) - 0.5f);
        } else {
            longitude = 0;
        }

        return Vector2f (latitude, longitude);
    }


    Point2f pixelPosition (const BoundingBox2i &dataWindow, const Point2f &latLong)
    {

        float x = latLong(1) / (-2 * M_PI) + 0.5f;
        float y = latLong(0) / -M_PI + 0.5f;


        return Vector2f (x * (dataWindow.max(0) - dataWindow.min(0)) + dataWindow.min(0),
            y * (dataWindow.max(1) - dataWindow.min(1)) + dataWindow.min(1));
    }


    Point2f pixelPosition (const BoundingBox2i &dataWindow, const Vector3f &direction)
    {
        return pixelPosition (dataWindow, latLong (direction));
    }


    Vector3f direction (const BoundingBox2i &dataWindow, const Point2f &pixelPosition){
        Vector2f ll = latLong (dataWindow, pixelPosition);

        return Vector3f (sin (ll(1)) * cos (ll(0)),
            sin (ll(0)),
            cos (ll(1)) * cos (ll(0)));
    }


};



class EnvMap : public Emitter {
public:
    EnvMap(const PropertyList &props) :
        Emitter(props.getTransform("toWorld", Transform()))
    {
        //set the arguments
        std::string filePath = props.getString("lightProbPath", "");

       filesystem::path path = getFileResolver()->resolve(filePath);
        //try to set it
        try {
            if (path.extension() == "exr") {
                m_lightprobe = new Bitmap(path.str());
                m_wrapper = new HSWrapper(*m_lightprobe);
                m_latlngMap = new LatLongMap();
                m_BBox.min = Point2i(0, 0);
                m_BBox.max = Point2i(m_lightprobe->cols(), m_lightprobe->rows());
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

    ~EnvMap(){
        delete m_wrapper;
        delete m_lightprobe;
        delete m_latlngMap;
    }

    std::string toString() const {
        return "envMap";
    }

    Color3f Power(const Scene *) const {return Color3f(0.0f); }

    bool isDeltaLight() const { return false; }

    Color3f sampleL(Point3f &p,
                    float pEpsilon,
                    const Point2f &ls, /*float time,*/
                    Vector3f *wi,
                    float *pdf,
                    VisibilityTester *vis) const {

        //set the direction to the light
        Point2f samplePix =  Warp::squareToHSW(ls, m_wrapper);
        Point2f normalizedPix  = samplePix;
        samplePix(0) *= m_width;
        samplePix(1) *= m_height;
        Vector3f sampleDirection = m_latlngMap->direction(m_BBox, samplePix);

        //std::cout << samplePix(0) << "; " << m_height - samplePix(1) << std::endl;

        //transform the sample to world coordinates
        *wi = LightToWorld * sampleDirection;        

        float sinTheta = sampleDirection(1);
        float imgPDF = Warp::squareToHSWPDF(normalizedPix, m_wrapper) ;
        //std::cout << imgPDF << std::endl;
        imgPDF *= m_height * m_width;
        //set the pdf of the sample
        *pdf =  (imgPDF ) / (sinTheta * 2.0f * M_PI * M_PI);

        //init the visibility tester
        vis->SetRay(p, *wi);


        Color3f result = m_wrapper->m_lightprob(samplePix(1), samplePix(0));
        return 0.0f < *pdf? result : Color3f(0.0);

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


    Color3f sampleL(const Vector3f &d) const {
        Vector3f localVector = WorldToLight * d;

        Point2f pixCoords = m_latlngMap->pixelPosition(m_BBox, localVector);

        int x = int(std::min(m_width  - 1.0f ,std::floor(pixCoords(0))));
        int y = int(std::min(m_height - 1.0f ,std::floor(pixCoords(1))));

        Color3f result(0.0f, 0.0f, 0.0f);
        bilinearInterpolated(x, y, result);

        return result;
    }

    Color3f radiance() const {
        return Color3f(0.0f);
    }
    void setMaxRadius(float radius){
        m_radius = radius;
    }

    float pdf(Point3f p, Vector3f w, Point3f hitPoint, Normal3f n) const {
        Vector3f localW = WorldToLight * w;
        Point2f pixCoord = m_latlngMap->pixelPosition(m_BBox, localW);
        return Warp::squareToHSWPDF(pixCoord, m_wrapper);
    }
    void samplePhoton(Sampler *sampler, Photon &photon, int N, int nLights, Vector3f &unQuantDir) const {
        Point2f samplePix =  Warp::squareToHSW(sampler->next2D(), m_wrapper);
        samplePix(0) *= m_width;
        samplePix(1) *= m_height;
        Vector3f sampleDirection = m_latlngMap->direction(m_BBox, samplePix);
        unQuantDir = sampleDirection;
        //cosine hemisphere normal
        Frame toHemiSphere(-sampleDirection);

        Vector3f cosDirection = Warp::squareToCosineHemisphere(sampler->next2D());

        Point3f pos = m_radius * toHemiSphere.toWorld(cosDirection);

        // hemisphre
        float A = 2.0f * M_PI * m_radius * m_radius ;

        Color3f power = A * m_wrapper->m_lightprob(samplePix(1), samplePix(0));

        // create the photon
        photon = Photon(pos, sampleDirection, power);
    }

private:
    HSWrapper *m_wrapper = nullptr;
    Bitmap *m_lightprobe = nullptr;
    LatLongMap *m_latlngMap = nullptr;
    BoundingBox2i m_BBox;
    int m_width = 0;
    int m_height = 0;
    float m_radius = 1.0f;

};

NORI_REGISTER_CLASS(EnvMap, "envMap");
NORI_NAMESPACE_END
