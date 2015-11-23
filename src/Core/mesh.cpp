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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/medium.h>
#include <nori/texture.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { }

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
    if(m_dpdf != nullptr) delete m_dpdf;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
    //check if the mesh is an emitter
    if(isEmitter()){
        m_dpdf = new DiscretePDF(m_F.cols());

        //for every face set the surface area
        for (int i = 0; i < m_F.cols(); ++i) {
            float curArea = surfaceArea(i);
            m_dpdf->append(curArea);
            m_surfaceArea += curArea;
        }

        //normalize the pdf
        m_dpdf->normalize();
    }
}

float Mesh::Pdf(const Point3f &p, const Point3f &hitP, const Normal3f &n, const Vector3f &wi) const {
    float distSquared = (hitP - p).squaredNorm();

    float AbsDot = std::abs(n.dot(- wi));


    if(AbsDot != 0.0f && m_surfaceArea != 0.0f) {
        return distSquared / (AbsDot * m_surfaceArea);
    } else {
        return 0.0f;
    }
}

void Mesh::samplePosition(const Point2f &sample, Point3f &p, Normal3f &n) const {
    Point2f reusableSample = sample;
    //get a random face
    uint32_t faceIndex = m_dpdf->sampleReuse(reusableSample(0));


    // sample the triangle of the face
    Point2f uv = Warp::squareToTriangle(reusableSample);

    //get the vertex index
    uint32_t i0 = m_F(0, faceIndex), i1 = m_F(1, faceIndex), i2 = m_F(2, faceIndex);

    //get the points
    Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    //compute the edge vectors
    Vector3f edge1 = (p1 - p0);
    Vector3f edge2 =  (p2 - p0);
    Vector3f e1 = uv(0) * edge1;
    Vector3f e2 = uv(1) * edge2;

    //compute the new point
    p = p0 + e1 + e2;

    //compute the normal base on baycentric coordinates
    //get the vertex normals
    if(m_N.cols() > 3) {
        Normal3f n0 = m_N.col(i0), n1 = m_N.col(i1), n2 = m_N.col(i2);
        Vector3f w0 = p0 - p;
        Vector3f w1 = p1 - p;
        Vector3f w2 = p2 - p;

        float area = (edge1.cross(edge2)).norm();
        if(area != 0.0f){
            float u = (w0.cross(w1)).norm() / area;
            float v = (w1.cross(w2)).norm() / area;

            n = (1.0f - u - v) * n0 + v * n1 + u * n2;
            n.normalized();

        } else {
            n = (e1.cross(e2)).normalized();

        }
    } else {
        //no normals provided
        n = (e1.cross(e2)).normalized();
        if(isnan(n.sum())){
            cout << "Normal is nan" << endl;
            cout << "e1 " << e1 << endl;
            cout << "e2 " << e2 << endl;
        }
    }




    //n = uv(0) * n0 + uv(1) * n1 + (1.0f - uv(0) - uv(1)) * n2;

}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;
        case ETEXTURE:
            if (m_texture)
                throw NoriException(
                    "Mesh: tried to register multiple texture instances!");
            m_texture = static_cast<Texture *>(obj);
        break;
        case EBUMPMAP:
            if (m_bumpmap)
                throw NoriException(
                    "Mesh: tried to register multiple bumpmap instances!");
            m_bumpmap = static_cast<BumpTexture *>(obj);
        break;
        case EMIXTEXTURE:
            if (m_texture)
                throw NoriException(
                    "Mesh: tried to register multiple texture instances!");
            m_texture = static_cast<Texture *>(obj);

        break;
        case EMIXBUMPMAP:
            if (m_bumpmap)
                throw NoriException(
                    "Mesh: tried to register multiple bumpmap instances!");
            m_bumpmap = static_cast<BumpTexture *>(obj);
        break;
        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;
        case EMedium: {

            Medium *media = static_cast<Medium *>(obj);
            if (m_medium)
                throw NoriException(
                    "Mesh: tried to register multiple medium instances!");
            m_medium = media;

        }
        break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
