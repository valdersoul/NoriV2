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

#include <nori/bitmap.h>
#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfVersion.h>
#include <ImfIO.h>
//#include <Magick++.h>
#include <png.hpp>

NORI_NAMESPACE_BEGIN

Bitmap::Bitmap(const std::string &filename) {
    Imf::InputFile file(filename.c_str());
    const Imf::Header &header = file.header();
    const Imf::ChannelList &channels = header.channels();

    Imath::Box2i dw = file.header().dataWindow();
    resize(dw.max.y - dw.min.y + 1, dw.max.x - dw.min.x + 1);

    cout << "Reading a " << cols() << "x" << rows() << " OpenEXR file from \""
         << filename << "\"" << endl;

    const char *ch_r = nullptr, *ch_g = nullptr, *ch_b = nullptr;
    for (Imf::ChannelList::ConstIterator it = channels.begin(); it != channels.end(); ++it) {
        std::string name = toLower(it.name());

        if (it.channel().xSampling != 1 || it.channel().ySampling != 1) {
            /* Sub-sampled layers are not supported */
            continue;
        }

        if (!ch_r && (name == "r" || name == "red" || 
                endsWith(name, ".r") || endsWith(name, ".red"))) {
            ch_r = it.name();
        } else if (!ch_g && (name == "g" || name == "green" || 
                endsWith(name, ".g") || endsWith(name, ".green"))) {
            ch_g = it.name();
        } else if (!ch_b && (name == "b" || name == "blue" || 
                endsWith(name, ".b") || endsWith(name, ".blue"))) {
            ch_b = it.name();
        }
    }

    if (!ch_r || !ch_g || !ch_b)
        throw NoriException("This is not a standard RGB OpenEXR file!");

    size_t compStride = sizeof(float),
           pixelStride = 3 * compStride,
           rowStride = pixelStride * cols();

    char *ptr = reinterpret_cast<char *>(data());

    Imf::FrameBuffer frameBuffer;
    frameBuffer.insert(ch_r, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert(ch_g, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert(ch_b, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); 
    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);

    m_totalLuminance = getTotalLuminace();
}

void Bitmap::save(const std::string &filename) {
    cout << "Writing a " << cols() << "x" << rows() 
         << " OpenEXR file to \"" << filename << "\"" << endl;

    Imf::Header header((int) cols(), (int) rows());
    header.insert("comments", Imf::StringAttribute("Generated by Nori"));

    Imf::ChannelList &channels = header.channels();
    channels.insert("R", Imf::Channel(Imf::FLOAT));
    channels.insert("G", Imf::Channel(Imf::FLOAT));
    channels.insert("B", Imf::Channel(Imf::FLOAT));

    Imf::FrameBuffer frameBuffer;
    size_t compStride = sizeof(float),
           pixelStride = 3 * compStride,
           rowStride = pixelStride * cols();

    char *ptr = reinterpret_cast<char *>(data());
    frameBuffer.insert("R", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert("G", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert("B", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); 

    Imf::OutputFile file(filename.c_str(), header);
    file.setFrameBuffer(frameBuffer);
    file.writePixels((int) rows());
}

void Bitmap::savePNG(const std::string &filename) {
    cout << "{\"message\":\"update\"}" << endl;
    png::image<png::rgb_pixel> img(cols(), rows());

        for(int i=0;i < rows();i++)
            for(int j=0;j<cols();j++)
            {
                Color4f  curColor = coeffRef(i, j);
                double r = std::pow(curColor(0), 1.0f / m_gamma);
                r = r > 1.0d ? 1.0d: r;
                r = r < 0.0d ? 0.0d: r;
                double g = std::pow(curColor(1), 1.0f / m_gamma);
                g = g > 1.0d ? 1.0d: g;
                g = g < 0.0d ? 0.0d: g;
                double b = std::pow(curColor(2), 1.0f / m_gamma);
                b = b > 1.0d ? 1.0d: b;
                b = b < 0.0d ? 0.0d: b;
                int iRed = (int) (r * 255.0d);
                int iGreen = (int) (g * 255.0d);
                int iBlue = (int) (b * 255.0d);
                png::rgb_pixel pixel(iRed, iGreen, iBlue);
                img.set_pixel(j, i, pixel);
            }


        //Write png to file
        img.write(filename);
    /*
    const Magick::Geometry size(cols(), rows());
    const Magick::ColorRGB color(0.0, 0.0, 0.0);
    Magick::Image img(size, color);

    // Allocate pixel view
    Magick::Pixels view(img);

    Magick::PixelPacket *pixels = view.get(0, 0, cols(), rows());

    for ( ssize_t row = 0; row < rows() ; ++row )
       for ( ssize_t column = 0; column < cols() ; ++column ) {
            Color4f  curColor = coeffRef(row, column);
            double r = std::pow(curColor(0), 1.0f / m_gamma);
            r = r > 1.0d ? 1.0d: r;
            r = r < 0.0d ? 0.0d: r;
            double g = std::pow(curColor(1), 1.0f / m_gamma);
            g = g > 1.0d ? 1.0d: g;
            g = g < 0.0d ? 0.0d: g;
            double b = std::pow(curColor(2), 1.0f / m_gamma);
            b = b > 1.0d ? 1.0d: b;
            b = b < 0.0d ? 0.0d: b;

            Magick::ColorRGB curColorMagick(r, g, b);
            *pixels++=curColorMagick;
       }


    // Save changes to image.
    view.sync();

    img.write(filename);
    */

}
float Bitmap::getTotalLuminace(){
    Color3f totalColor = sum();

    float totalLum = totalColor.getLuminance();
    return totalLum;
}


NORI_NAMESPACE_END
