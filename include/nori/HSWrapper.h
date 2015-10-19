#include <nori/bitmap.h>
#include <nori/sampler.h>
#include <nori/dpdf.h>
NORI_NAMESPACE_BEGIN

struct imagQuad {
    Point2i BBMin;
    Point2i BBMax;
    float propability = 0.0f;
    std::unique_ptr<imagQuad> child1;
    std::unique_ptr<imagQuad> child2;
    std::unique_ptr<imagQuad> child3;
    std::unique_ptr<imagQuad> child4;
    bool leave = false;

    imagQuad() { }

    imagQuad(Point2i start, Point2i stop, Bitmap &lightprob) {
        // set the bounding box
        BBMin = start;
        BBMax = stop;
        //this quad part is a leave
        if(BBMax(0) - BBMin(0) < 2 || BBMax(1) - BBMin(1) < 2) {
            leave = true;
            return;
        } else {
            //find the middle
            Point2i middle = BBMin;
            middle(0) += int(std::floor((BBMax(0) - BBMin(0)) / 2.0));
            middle(1) += int(std::floor((BBMax(1) - BBMin(1)) / 2.0));

            //get the luminance of the image
            /*
             * __________________
             * |        |        |
             * |   L1   |  L2    |
             * |        |        |
             * |--------|--------|
             * |   L3   |  L4    |
             * |        |        |
             * |________|________|
             */
            int i = BBMin(0);
            int j = BBMin(1);
            int l = middle(0);
            int k = middle(1);
            Bitmap lightProb1(lightprob.block(i, j, l, k));
            i = BBMin(0);
            j = middle(1);
            l = middle(1) - i;
            k = BBMax(1) - j;
            Bitmap lightProb2(lightprob.block(i, j, l, k));
            i = middle(0);
            j = BBMin(1);
            l = BBMax(0) - i;
            k = middle(1) - j;
            Bitmap lightProb3(lightprob.block(i, j, l, k));
            i = middle(0);
            j = middle(1);
            l = BBMax(0) - i;
            k = BBMax(1) - j;
            Bitmap lightProb4(lightprob.block(i, j, l, k));

            //calculate the probabilities
            float totalLuminace1 = lightProb1.getTotalLuminace();
            float totalLuminace2 = lightProb2.getTotalLuminace();
            float totalLuminace3 = lightProb3.getTotalLuminace();
            float totalLuminace4 = lightProb4.getTotalLuminace();
            float totalLuminace = totalLuminace1 + totalLuminace2 + totalLuminace3 + totalLuminace4;

            //init the childeren
            child1.reset(new imagQuad(Point2i(BBMin(0), BBMin(1)), Point2i(middle(0), middle(1)), lightProb1));
            child2.reset(new imagQuad(Point2i(BBMin(0), middle(1)), Point2i(middle(1), BBMax(1)), lightProb2));
            child3.reset(new imagQuad(Point2i(middle(0), BBMin(1)), Point2i(BBMax(0), middle(1)), lightProb3));
            child4.reset(new imagQuad(Point2i(middle(0), middle(1)), Point2i(BBMax(0), BBMax(1)), lightProb4));

            // set the probalitiy
            if(totalLuminace != 0.0f) {
                child1->propability = totalLuminace1 / totalLuminace;
                child2->propability = totalLuminace2 / totalLuminace;
                child3->propability = totalLuminace3 / totalLuminace;
                child4->propability = totalLuminace4 / totalLuminace;
            }
        }
    }
    void sample(Point2i &pixIndex, Point2f sample) {
        if(leave){
            //get the pixel postion
            pixIndex(0) = int(std::floor((BBMax(0) - BBMin(0)) * sample.y()));
            pixIndex(1) = int(std::floor((BBMax(1) - BBMin(1)) * sample.x()));
            return;
        }
        //check if upper lower childern
        float probY = child1->propability + child2->propability;
        if(sample.y() < probY){
            //check upper left or right
            //resample y
            sample.y() = sample.x() / probY;
            if(sample.x() < child1->propability) {
                //resample x
                sample.x() = sample.x() / child1->propability;
                child1->sample(pixIndex, sample);
            } else {
                sample.x() = (sample.x() - child1->propability) / (1.0f - child1->propability);
                child2->sample(pixIndex, sample);
            }
        } else {
            //check lower left or right
            sample.y() = (sample.y() - probY) / (1.0f - probY);
            if(sample.x() < child3->propability) {
                sample.x() = sample.x() / child2->propability;
                child3->sample(pixIndex, sample);
            } else {
                sample.x() = (sample.x() - child3->propability) / (1.0f - child3->propability);
                child4->sample(pixIndex, sample);
            }
        }

    }

};
class HSWrapper {
public:
    HSWrapper(){}
    ~HSWrapper(){
        if(!m_useHSW){
            delete m_dpdfRow;
            for (uint32_t i = 0; i < colDPdfs.size(); ++i) {
                delete colDPdfs[i];
            }
        } else {
            delete m_root;
        }

    }

    HSWrapper(Bitmap &curLightprob){
        //set the current light prob
        m_lightprob = curLightprob;

        if(m_useHSW){
            m_root = new imagQuad(Point2i(0, 0), Point2i(curLightprob.rows(), curLightprob.cols()), curLightprob);
        } else {
            m_dpdfRow = new DiscretePDF(m_lightprob.rows());

            //for every row the image
            for (uint32_t i = 0; i < m_lightprob.rows(); ++i) {
                float totalRowLum = 0.0f;
                //add up all pixels in the row
                DiscretePDF *curColPdf = new DiscretePDF(m_lightprob.cols());
                for (uint32_t j = 0; j < m_lightprob.cols(); ++j) {
                    Color3f curColor = m_lightprob(i, j);
                    float curLum = curColor.getLuminance();
                    curColPdf->append(curLum);
                    totalRowLum += curLum;
                }

                curColPdf->normalize();
                colDPdfs.push_back(curColPdf);

                m_dpdfRow->append(totalRowLum);
            }

            //normalize the pdf
            m_dpdfRow->normalize();
        }


        for (int y = 0; y < m_lightprob.rows(); ++y) {
            for (int x = 0; x < m_lightprob.cols(); ++x) {
                m_totalLum  += m_lightprob(y, x).getLuminance();
            }

        }

    }

    void sample(Point2f sample, Color3f &result, Point2i &pixCoords);

    Bitmap m_lightprob;
    float m_totalLum = 0.0f;
    DiscretePDF *m_dpdfRow = nullptr;
    std::vector<DiscretePDF*> colDPdfs;
private:
    imagQuad *m_root = nullptr;
    bool m_useHSW = false;


};

NORI_NAMESPACE_END
