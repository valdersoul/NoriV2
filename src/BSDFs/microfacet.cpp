#include <nori/microfacet.h>

NORI_NAMESPACE_BEGIN

float Microfacet::G1(float alpha, const Vector3f &v, const Vector3f &m){
    if (v.dot(m) * Frame::cosTheta(v) <= 0.0f)
        return 0.0f;

    float c = (v.dot(m)) / (Frame::cosTheta(v));
    float b = 1.0f / (alpha * Frame::tanTheta(v));

    if(b < 1.6) {
        return CHIPlus(c) * ((3.535f * b + 2.181f * b * b) / (1.0 + 2.276 * b + 2.577 * b * b));
    }

    return CHIPlus(c);
}
float Microfacet::CHIPlus(const float c) {
    return 0.0f < c ? 1.0f : 0.0f;
}

float Microfacet::D(float alpha, const Vector3f &m){
    if (m.z() <= 0.0f)
        return 0.0f;

    return Warp::squareToBeckmannPdf(m, alpha);
}


float Microfacet::G(float alpha, const Vector3f &i, const Vector3f &o, const Vector3f &m){
    return G1(alpha, i, m) * G1(alpha, o, m);
}

float Microfacet::pdf(float alpha, const Vector3f &m){
    return Microfacet::D(alpha, m) * Frame::cosTheta(m);
}

Vector3f Microfacet::sample(float alpha, Point2f xi){
    return Warp::squareToBeckmann(xi, alpha);
}

NORI_NAMESPACE_END
