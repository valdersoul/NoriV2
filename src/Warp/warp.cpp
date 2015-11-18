/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <algorithm>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    //get the radius (uniform)
    float r = sqrt(sample(0));
    //get the theta angle
    float theta = 2.0f * M_PI *  sample(1);

    //create a point on the disk
    Point2f res(r * cos(theta), r * sin(theta));
    //cout << "Sampled point = " << res << endl;
    return res;
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    //check if point is inside the unitdisk
    return (p.norm() < 1) ? INV_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    //convert cos theta max to a hight of the cylinder
    float hightStop= (1.0 - cosThetaMax) / 2.0;

    Point2f cylSample = squareToUniformCylinder(sample, cosThetaMax, hightStop);

    float omegaZ = cylSample(0);
    float r = sqrt(1.0 -  omegaZ * omegaZ);
    float phi = cylSample(1);
    float omegaX = r * cos(phi);
    float omegaY = r * sin(phi);

    return Vector3f(omegaX, omegaY, omegaZ);
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {    
    if(cosThetaMax != 1.0) {
        return (v(2) < cosThetaMax) ? 0.0 : 1.0 / (2.0 * M_PI * (1.0f - cosThetaMax));
    } else {
        if(v(2) == 1.0) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
}

Point2f Warp::squareToHSW(const Point2f &sample, HSWrapper *testWrapper) {
    Color3f result;
    Point2i pixCoord;
    testWrapper->sample(sample, result, pixCoord);
    int height = testWrapper->m_lightprob.rows();
    int width = testWrapper->m_lightprob.cols();
    float y = float(pixCoord(0)) / float(height);
    float x = float(pixCoord(1)) / float(width);
    Point2f pixResult(x, y);

    return pixResult;
}
float Warp::squareToHSWPDF(const Point2f &p, HSWrapper *testWrapper) {

    float total = testWrapper->m_totalLum;
    float height = testWrapper->m_lightprob.rows();
    float width = testWrapper->m_lightprob.cols();
    float fPixH = height * p(1);
    float fPixW = width * p(0);
    int y = int(std::min(height - 1.0f ,std::floor(fPixH)));
    int x = int(std::min(width  - 1.0f ,std::floor(fPixW)));

    int yNext = y + 1;
    int xNext = x + 1;

    Color3f Q11 = testWrapper->m_lightprob(y, x);
    Color3f Q21 = testWrapper->m_lightprob(y, std::min(int(width -1.0f), xNext));
    Color3f Q12 = testWrapper->m_lightprob(std::min(int(height -1.0f), yNext), x);
    Color3f Q22 = testWrapper->m_lightprob(std::min(int(height -1.0f), yNext), std::min(int(width -1.0f), xNext));

    float fact1 = ((xNext - fPixW) / (xNext - x));
    float fact2 = ((fPixW - x) / (xNext - x));

    Color3f R1 = fact1 * Q11 + fact2 * Q21;
    Color3f R2 = fact1 * Q12 + fact2 * Q22;

    fact1 = ((yNext - fPixH) / (yNext - y));
    fact2 = ((fPixH - y) / (yNext - y));

    Color3f result = fact1 * R1 + fact2 * R2;

    float lum = result.getLuminance();
    if (lum < 0) {
        cout << "Negative luminance detected : " << lum << endl;
    }
    float pdf = lum / total;

    return pdf;

}
Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    Point2f cylSample = squareToUniformCylinder(sample, -1.0f, 1.0f);
    float omegaZ = cylSample(0);
    float r = sqrt(1.0 -  omegaZ * omegaZ);
    float phi = cylSample(1);
    float omegaX = r * cos(phi);
    float omegaY = r * sin(phi);
    return Vector3f(omegaX, omegaY, omegaZ);


}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return (v.norm() <= 1.0) ? INV_FOURPI : 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    Point2f cylSample = Warp::squareToUniformCylinder(sample, 0.0, 0.5);
    float omegaZ = cylSample(0);
    float r = sqrt(1.0 -  omegaZ * omegaZ);
    float phi = cylSample(1);
    float omegaX = r * cos(phi);
    float omegaY = r * sin(phi);
    return Vector3f(omegaX, omegaY, omegaZ);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    //return (v.norm() <= 1.0 && v(2) >= 0.0) ? INV_TWOPI : 0.0f;
    return (v(2) >= 0.0) ? 1.0 / (2.0 * M_PI) : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    //create a sample which is cosine weighted
   Point2f ptOnDisk =  squareToUniformDisk(sample);

   //compute the Z value
   return Vector3f(ptOnDisk(0), ptOnDisk(1), sqrt(1.0 - ptOnDisk.squaredNorm()));


}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return (v(2) >= 0) ? v(2) * INV_PI : 0.0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {

    //use the inverse methode

    float theta = sqrt(-1.0 / (alpha * alpha * log(sample(0)) -1.0 ));
    float phi = 2.0 * M_PI * sample(1);

    theta = acos(theta);
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

}


float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    //compute the angle
    float theta = acos(m(2));
    if(theta < M_PI_2) {
        //float phi = atan(m(1) / m(0));

        //calculate the beckmann pdf
        float cosTheta2 = cos(theta) * cos(theta);
        float bEXP = (m.x() * m.x()) / (alpha * alpha)
                + (m.y()*m.y() / (alpha * alpha)) / (cosTheta2);

        float num = exp(- bEXP);
        float denum = M_PI * alpha * alpha * cosTheta2 * cosTheta2;

        return num / denum ;
    }
    return 0.0;
}

Point2f Warp::squareToUniformCylinder(const Point2f &sample, float hightStart, float hightStop) {
    return Point2f(2.0f * hightStop * sample(1) + hightStart, 2.0F * M_PI  * sample(0));
}

Point2f Warp::squareToTriangle(const Point2f &sample){
    float s = 1.0f - sqrt(1.0f - sample(0));
    float t = (1.0f - s) * sample(1);
    return Point2f(std::max(s, 0.0001f), std::max(t, 0.0001f));
}

Vector3f Warp::squareToHenyeyGreenstein(const Point2f &sample, float g) {

    //use the inverse methode
    if(g == 0.0){
        return Warp::squareToUniformSphere(sample);
    }
    float scale = 1.0f / (2.0f * g);
    float fraction = (1.0f - g*g) / (1.0f - g + 2.0f * g * std::max(sample(0), 0.001f));
    float theta =  scale * (1.0f +  g * g - fraction * fraction);


    float phi = 2.0 * M_PI * sample(1);

    float acostheta = acos(theta);

    return Vector3f(sin(acostheta) * cos(phi), sin(acostheta) * sin(phi), cos(acostheta));
}
float Warp::squareToHenyeyGreensteinPDF(const Vector3f &p, float g){
    float cosTheta = p(2);
    float sqrtFrac = (1.0f + g * g - 2.0f * g * cosTheta);
    if(sqrtFrac <= 0.0f)
        return 0.0f;
    float pdf =  ((1.0f - g * g) / (4.0f * M_PI * pow(sqrtFrac, 1.5f)));

    return pdf;
}
float Warp::squareToHenyeyGreensteinPDF(const Vector3f &wo, const Vector3f &wi, float g){
    float cosTheta = wo.dot(wi);
    float sqrtFrac = (1.0f + g * g - 2.0f * g * cosTheta);
    if(sqrtFrac <= 0.0f)
        return 0.0f;
    float pdf =  ((1.0f - g * g) / (4.0f * M_PI * pow(sqrtFrac, 1.5f)));

    return pdf;
}

Vector3f Warp::squareToSchlick(const Point2f &sample, float k) {

    //use the inverse methode
    if(k == 0.0){
        return Warp::squareToUniformSphere(sample);
    }

    float cosTheta = ((k + 2.0f * sample.x() - 1.0f ) / ( k * (2.0f * sample.x() - 1.0f ) + 1.0f));
    float acostheta = acos(cosTheta);
    float phi = 2.0f * M_PI * sample.y();

    return Vector3f(sin(acostheta) * cos(phi), sin(acostheta) * sin(phi), cos(acostheta));
}

float Warp::squareToSchlickPDF(const Vector3f &p, float k){
    float cosTheta = p(2);
    float num = 1.0f - k * k;
    float denum = 4.0f * M_PI * (1.0f -  k * cosTheta) * (1.0f -  k * cosTheta);
    float pdf = num / denum;
    return pdf;
}

float Warp::squareToSchlickPDF(const Vector3f &wo, const Vector3f &wi, float k){
    float cosTheta = wo.dot(wi);
    float num = 1.0f - k * k;
    float denum = 4.0f * M_PI * (1.0f -  k * cosTheta) * (1.0f -  k * cosTheta);
    float pdf = num / denum;
    return pdf;
}

NORI_NAMESPACE_END
