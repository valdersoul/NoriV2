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

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <nori/areaLight.h>
//#include <nori/volume.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &props) {
    m_bvh = new BVH();
    m_numberOfPasses = props.getInteger("nrPasses", 1);
}

Scene::~Scene() {
    delete m_bvh;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
}

void Scene::activate() {
    m_bvh->build();

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

void Scene::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_bvh->addMesh(mesh);
                m_meshes.push_back(mesh);

                //check if the mesh is an emitter
                if(mesh->isEmitter()) {
                    Emitter* areaLightEM = mesh->getEmitter();

                    //add the Emitter to the list of emitters
                    m_emitters.push_back(areaLightEM);

                    areaLight* aEM = static_cast<areaLight *>(areaLightEM);
                    aEM->setMesh(mesh);
                }
            }
            break;
        
        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                m_emitters.push_back(emitter);

                if(emitter->toString() == "diskLight" || emitter->toString() == "envMap") {
                    m_distantEmitter = emitter;
                }

                break;
                //throw NoriException("Scene::addChild(): You need to implement this for emitters");
            }
            break;

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;
        case EVolume:
            {
                Volume *vol = static_cast<Volume *>(obj);
                m_volumes.push_back(vol);
                break;
            }
            break;
        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    std::string volumes;
    for (size_t i=0; i<m_volumes.size(); ++i) {
        volumes += std::string("  ") + indent(m_volumes[i]->toString(), 2);
        if (i + 1 < m_volumes.size())
            volumes += ",";
        volumes += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  volumes = {\n"
        "  %s  }\n"
        "  meshes = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(volumes, 2),
        indent(meshes, 2)
    );
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
