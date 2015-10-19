#include <nori/HSWrapper.h>

NORI_NAMESPACE_BEGIN


void HSWrapper::sample(Point2f sample, Color3f &result, Point2i &pixCoords){

    if(m_useHSW) {
        // get postion in the image
        m_root->sample(pixCoords, sample);
    } else {

        //DiscretePDF dpdfCol(m_lightprob.cols());

        //sample a row
        int rowIndex = m_dpdfRow->sample(sample(0));

        //sample the column
        //int colIndex = dpdfCol.sample(sample(1));
        int colIndex = colDPdfs[rowIndex]->sample(sample(1));

        //the the resulting pixel coordiantes
        pixCoords(0) = rowIndex;
        pixCoords(1) = colIndex;
    }

    //get the color of the envmap
    result = m_lightprob(pixCoords(0), pixCoords(1));
}

NORI_NAMESPACE_END
