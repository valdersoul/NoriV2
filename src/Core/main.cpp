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

#include <nori/parser.h>
#include <nori/scene.h>
#include <nori/camera.h>
#include <nori/block.h>
#include <nori/timer.h>
#include <nori/bitmap.h>
#include <nori/sampler.h>
#include <nori/integrator.h>
#include <nori/gui.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>
#include <filesystem/resolver.h>
#include <thread>
#include <tbb/task_scheduler_init.h>
#include <stdlib.h>

using namespace nori;
typedef tbb::queuing_mutex Mutex;


static void renderBlock(const Scene *scene, Sampler *sampler, ImageBlock &block) {
    const Camera *camera = scene->getCamera();
    const Integrator *integrator = scene->getIntegrator();

    Point2i offset = block.getOffset();
    Vector2i size  = block.getSize();

    /* Clear the block contents */
    block.clear();

    /* For each pixel and pixel sample sample */
    for (int y=0; y<size.y(); ++y) {
        for (int x=0; x<size.x(); ++x) {
            for (uint32_t i=0; i<sampler->getSampleCount(); ++i) {
                Point2f pixelSample = Point2f((float) (x + offset.x()), (float) (y + offset.y())) + sampler->next2D();
                Point2f apertureSample = sampler->next2D();

                /* Sample a ray from the camera */
                Ray3f ray;
                Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);

                /* Compute the incident radiance */
                value *= integrator->Li(scene, sampler, ray);

                /* Store in the image block */
                block.put(pixelSample, value);
            }
        }
    }
}

static void render(Scene *scene, const std::string &filename, const int nrThreads, const bool showWindow, const bool saveEveryStep) {
    const Camera *camera = scene->getCamera();
    Vector2i outputSize = camera->getOutputSize();

    NoriScreen *screen = nullptr;

    /* Allocate memory for the entire output image and clear it */
    ImageBlock result(outputSize, camera->getReconstructionFilter());
    result.clear();

    /* Create a window that visualizes the partially rendered result */
    if (showWindow) {
        nanogui::init();
        screen = new NoriScreen(result);
    }

    /* Determine the filename of the output bitmap */
    std::string outputName = filename;
    std::string outputNamePNG = filename;
    size_t lastdot = outputName.find_last_of(".");
    if (lastdot != std::string::npos) {
        outputName.erase(lastdot, std::string::npos);
        outputNamePNG.erase(lastdot, std::string::npos);
    }
    outputName += ".exr";
    outputNamePNG += ".png";

    /* Do the following in parallel and asynchronously */


    std::thread render_thread([&] {
        cout << "Rendering .. ";
        cout.flush();
        Timer timer;

        do {
            scene->getIntegrator()->preprocess(scene);

            int nPasses = scene->getNumberOfPasses();

            //BEGIN MULTIPLE PASS LOOP - PROGRESSIVE RENDERING WITH THE SAME INTEGRATOR SETTINGS
            for (int p = 1; p <= nPasses; p++){
                std::cout << "Rendering pass " << p << "/" << nPasses << std::endl;
                /* Create a block generator (i.e. a work scheduler) */
                BlockGenerator blockGenerator(outputSize, NORI_BLOCK_SIZE);

                tbb::task_scheduler_init init(nrThreads);
                int blockCount = blockGenerator.getBlockCount();
                int percent = 0;
                tbb::blocked_range<int> range(0, blockCount);

                auto map = [&](const tbb::blocked_range<int> &range) {
                    /* Allocate memory for a small image block to be rendered
                       by the current thread */
                    ImageBlock block(Vector2i(NORI_BLOCK_SIZE),
                        camera->getReconstructionFilter());
                    Mutex mutex;



                    /* Create a clone of the sampler for the current thread */
                    std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());

                    for (int i=range.begin(); i<range.end(); ++i) {
                        /* Request an image block from the block generator */
                        blockGenerator.next(block);

                        /* Inform the sampler about the block to be rendered */
                        sampler->prepare(block, p);

                        /* Render all contained pixels */
                        renderBlock(scene, sampler.get(), block);

                        /* The image block has been processed. Now add it to
                           the "big" block that represents the entire image */
                        result.put(block);
                        if(saveEveryStep) {
                            Mutex::scoped_lock lock(mutex);

                            /* Save using the png format */
                            std::unique_ptr<Bitmap> bitmap(result.toBitmap());
                            bitmap->savePNG(outputNamePNG, int(100.0f * float(percent) / float(blockCount)));

                            std::unique_ptr<Bitmap> blockBitmap(block.toBitmap());

                            /* print the data as base64 string*/
                            blockBitmap->printB64Encoded(result.getSize(), block.getSize(), block.getOffset(), int(100.0f * float(percent) / float(blockCount)));

                            percent++;
                        }
                    }
                };

                /// Uncomment the following line for single threaded rendering
                // map(range);

                /// Default: parallel rendering
                tbb::parallel_for(range, map);
                /* Now turn the rendered image block into
                   a properly normalized bitmap */
                std::unique_ptr<Bitmap> bitmap(result.toBitmap());

                /* Save using the OpenEXR format */
                bitmap->save(outputName);
            }
        // DOES THE RENDERER WANT TO CONTINUE WITH DIFFERENT SETTINGS?
        // EXAMPLE - RADIUS CHANGE FOR PHOTON MAPPING
        } while (scene->getIntegrator()->advance());

        cout << "done. (took " << timer.elapsedString() << ")" << endl;
    });
    if(showWindow) {
        /* Enter the application main loop */
        nanogui::mainloop();
    }

    /* Shut down the user interface */
    render_thread.join();

    if(showWindow) {
        delete screen;
        nanogui::shutdown();
    }
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 5) {
        cerr << "Syntax: " << argv[0] << " <scene.xml> [number of threads] [show window (0|1)] [save every image (0|1)]" << endl;
        return -1;
    }

    filesystem::path path(argv[1]);

    try {
        if (path.extension() == "xml") {
            /* Add the parent directory of the scene file to the
               file resolver. That way, the XML file can reference
               resources (OBJ files, textures) using relative paths */
            getFileResolver()->prepend(path.parent_path());

            std::unique_ptr<NoriObject> root(loadFromXML(argv[1]));
            int numberOfThreads = 0;
            bool showWindow = true;
            bool saveEveryImage = false;
            char *endptr;
            if(argc > 2){
                numberOfThreads = int(std::strtol(argv[2], &endptr, 10));

            }
            if (argc > 3) {
                if(int(std::strtol(argv[3], &endptr, 10)) == 0) {
                    showWindow = false;
                }
            }
            if (argc > 4) {
                if(int(std::strtol(argv[4], &endptr, 10)) == 1) {
                    saveEveryImage = true;
                }
            }
            /* When the XML root object is a scene, start rendering it .. */
            if (root->getClassType() == NoriObject::EScene) {
                render(static_cast<Scene *>(root.get()), argv[1], numberOfThreads, showWindow, saveEveryImage);
            }
        } else if (path.extension() == "exr") {
            /* Alternatively, provide a basic OpenEXR image viewer */
            Bitmap bitmap(argv[1]);
            ImageBlock block(Vector2i(bitmap.cols(), bitmap.rows()), nullptr);
            block.fromBitmap(bitmap);
            nanogui::init();
            NoriScreen *screen = new NoriScreen(block);
            nanogui::mainloop();
            delete screen;
            nanogui::shutdown();
        } else {
            cerr << "Fatal error: unknown file \"" << argv[1]
                 << "\", expected an extension of type .xml or .exr" << endl;
        }
    } catch (const std::exception &e) {
        cerr << "Fatal error: " << e.what() << endl;
        return -1;
    }
    return 0;
}
