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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

#include <iostream>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the
           specular component by 1-kd.

           While that is not a particularly realistic model of what
           happens in reality, this will greatly simplify the
           implementation. Please see the course staff if you're
           interested in implementing a more realistic version
           of this BRDF. */
        m_ks = 1.0f - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {

        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        float cosThetai = Frame::cosTheta(bRec.wi);
        float cosThetao = Frame::cosTheta(bRec.wo);


        //compute the Beckman term
        float D = Warp::squareToBeckmannPdf(wh, m_alpha) / Frame::cosTheta(wh);

        //compute the Fresnel term
        float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);

        //compute the geometry term
        float G = G1(bRec.wi, wh) * G1(bRec.wo, wh);

        return (m_kd * INV_PI) + (m_ks * ((D * F * G)/(4.0f * cosThetai * cosThetao)));
    }

    float G1(Vector3f w, Vector3f wh) const{
        float c = (w.dot(wh)) / (Frame::cosTheta(w));
        //float b = 1.0f / (m_alpha * tan(n.dot(w)));
        float b = 1.0f / (m_alpha * Frame::tanTheta(w));

        if(b < 1.6) {
            return CHIPlus(c) * ((3.535 * b + 2.181 * b * b) / (1.0 + 2.276 * b + 2.577 * b * b));
        }
        return CHIPlus(c);
    }

    float CHIPlus(float c) const {
        return 0.0f < c ? 1.0f : 0.0f;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        //compute the Beckman term devided by another cosine
        float D = Warp::squareToBeckmannPdf(wh, m_alpha) / Frame::cosTheta(wh);

        float cosThetah = Frame::cosTheta(wh);
        float cosThetao = Frame::cosTheta(bRec.wo);

        //Jacobian of the half direction mapping
        float J = 1.0f / (4.0f * (wh.dot(bRec.wo)));

        float term1 = m_ks * D * cosThetah * J;
        float term2 = (1.0f - m_ks) * cosThetao * INV_PI;

        //cout << "term1 = " << term1 << endl;
        //  cout << "term2 = " << term2 << endl;

        return term1 + term2;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);
        bRec.measure = ESolidAngle;
        /* Warp a uniformly distributed sample on [0,1]^2
        to a direction on a cosine-weighted hemisphere */
        if (sample.x() < m_ks) {
            float x = sample.x()/m_ks;
            Vector3f n = Warp::squareToBeckmann(Point2f(x, sample.y()), m_alpha);
            bRec.wo = (2.0f * bRec.wi.dot(n) * n - bRec.wi );

        }
        else {
            float x = (sample.x() - m_ks) / (1.0f - m_ks);
            bRec.wo = Warp::squareToCosineHemisphere(Point2f(x,sample.y()));
        }
        if(Frame::cosTheta(bRec.wo) <= 0) {
            return Color3f(0.0f);
        }
        //std::cout << "sample.x = " << sample.x() << std::endl;
        //std::cout << "wo = " << bRec.wo.toString() << std::endl;

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;
        float pdfVal = pdf(bRec);
        return (pdfVal != 0.0f) ? Color3f(eval(bRec) * Frame::cosTheta(bRec.wo) / pdfVal) : Color3f(0.0f);

    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
