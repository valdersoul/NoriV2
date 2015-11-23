# This is the main CMake build file used to compile Nori
cmake_minimum_required (VERSION 2.8)
project("nori")

# This file included below takes care of compiling dependencies
# and setting compiler flags on different platforms
include("CMakeConfig.txt")

# The following lines build the main executable. If you add a source
# code file to Nori, be sure to include it in this list.
add_executable(nori

  # Header files
  include/nori/bumpTexture.h
  include/nori/mixBumpmap.h
  include/nori/mixTexture.h
  include/nori/perlinTexture.h
  include/nori/texture.h
  include/nori/HSWrapper.h
  include/nori/phaseFunction.h
  include/nori/volume.h
  include/nori/metalConsts.h

  include/nori/areaLight.h
  include/nori/bbox.h
  include/nori/bitmap.h
  include/nori/block.h
  include/nori/bumpTexture.h
  include/nori/bsdf.h
  include/nori/medium.h
  include/nori/bvh.h
  include/nori/camera.h
  include/nori/color.h
  include/nori/common.h
  include/nori/dpdf.h
  include/nori/frame.h
  include/nori/integrator.h
  include/nori/emitter.h
  include/nori/mesh.h
  include/nori/mixTexture.h
  include/nori/mixBumpmap.h
  include/nori/object.h
  include/nori/parser.h
  include/nori/proplist.h
  include/nori/ray.h
  include/nori/rfilter.h
  include/nori/sampler.h
  include/nori/scene.h
  include/nori/perlinTexture.h
  include/nori/timer.h
  include/nori/texture.h
  include/nori/transform.h
  include/nori/vector.h
  include/nori/warp.h
  include/nori/microfacet.h

  # Source code files

  # BSDFS
  src/BSDFs/diffuse.cpp
  src/BSDFs/microfacet.cpp
  src/BSDFs/microfacetBRDF.cpp
  src/BSDFs/mirror.cpp
  src/BSDFs/dielectric.cpp
  src/BSDFs/conductor.cpp
  src/BSDFs/roughconductor.cpp
  src/BSDFs/roughdielectric.cpp

  #medium
  src/Medium/homogeneous.cpp

  # Cameras
  src/Cameras/perspective.cpp

  # Core
  src/Core/bitmap.cpp
  src/Core/block.cpp
  src/Core/bvh.cpp
  src/Core/chi2test.cpp
  src/Core/common.cpp
  src/Core/gui.cpp
  src/Core/main.cpp
  src/Core/mesh.cpp
  src/Core/obj.cpp
  src/Core/object.cpp
  src/Core/parser.cpp
  src/Core/proplist.cpp
  src/Core/rfilter.cpp
  src/Core/scene.cpp
  src/Core/ttest.cpp
  src/Core/photon.cpp

  # Emitters
  src/Emitters/areaLight.cpp
  src/Emitters/pointLight.cpp
  src/Emitters/envMapEmitter.cpp
  src/Emitters/diskLight.cpp
  src/Emitters/spotLight.cpp
  src/Emitters/textureEmitter.cpp

  # Integrators
  src/Integrators/avgVisibility.cpp
  src/Integrators/emissionIntegrator.cpp
  src/Integrators/testtexture.cpp
  src/Integrators/testtexturePath.cpp
  src/Integrators/basicVolumentricPathTracer.cpp
  src/Integrators/direct_ems.cpp
  src/Integrators/direct_mats.cpp
  src/Integrators/direct_mis.cpp
  src/Integrators/direct.cpp
  src/Integrators/normals.cpp
  src/Integrators/path_mats.cpp
  src/Integrators/path_mis.cpp
  src/Integrators/photonmapper.cpp
  src/Integrators/testtexturePhoton.cpp
  src/Integrators/VPTWithShadowRays.cpp
  src/Integrators/volumetricPhotonmapper.cpp
  src/Integrators/beamRadianceEstimate.cpp
  src/Integrators/beamRESurfKNN.cpp
  src/Integrators/basicMediumPathTracer.cpp




  # Sampler
  src/Sampler/independent.cpp

  # Texture
  src/Texture/singleBumpTexture.cpp
  src/Texture/perlinTexture.cpp
  src/Texture/imageTexture.cpp
  src/Texture/mixTexture.cpp
  src/Texture/mixBumpmap.cpp

  # Volumes
  src/Volumes/homogeneousVolumeDensity.cpp

  # Warp
  src/Warp/HenyeyGreenstein.cpp
  src/Warp/HSWrapper.cpp
  src/Warp/warp.cpp
  src/Warp/Schlick.cpp
  

)

# The following lines build the warping test application
add_executable(warptest
  include/nori/warp.h
  include/nori/bitmap.h
  include/nori/HSWrapper.h

  src/Warp/HSWrapper.cpp
  src/Core/bitmap.cpp
  src/Warp/warp.cpp
  src/Warp/warptest.cpp
  src/BSDFs/microfacet.cpp
  src/Core/object.cpp
  src/Core/proplist.cpp
  src/Core/common.cpp
)

# Nori depends on some libraries created in CMakeConfig.txt. The following two
# lines ensure that Nori is built *after* those libraries have been created.
add_dependencies(nori OpenEXR_p)
add_dependencies(nori nanogui_p)
add_dependencies(nori tbb_p)
add_dependencies(nori pugixml)
add_dependencies(warptest nori)

# Link to dependency libraries
target_link_libraries(nori ${extra_libs})
target_link_libraries(warptest ${extra_libs})

# find image magick++
find_package(ImageMagick COMPONENTS Magick++)
if (ImageMagick_FOUND)
    include_directories(${ImageMagick_INCLUDE_DIRS})
    target_link_libraries (nori ${ImageMagick_LIBRARIES})
    target_link_libraries (warptest ${ImageMagick_LIBRARIES})
    MESSAGE( STATUS "ImageMagick found!")
endif (ImageMagick_FOUND)

# vim: set et ts=2 sw=2 ft=cmake nospell:
