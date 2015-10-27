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

#if !defined(__NORI_METALCONSTS_H)
#define __NORI_METALCONSTS_H

#include <nori/color.h>

NORI_NAMESPACE_BEGIN
/**
 * \brief Convenience data structure used to store multiple
 * parameters for conductors
 */
struct metalMaterial {
    // the Chemical material name
    std::string name;
    // the refraction coefficient
    Color3f eta;
    // the absorbtion coefficient
    Color3f k;

    // the constructor
    metalMaterial(const std::string name, const Color3f eta, const Color3f k)
        : name(name), eta(eta), k(k) { }
};
static const metalMaterial matrialList[] = {
            {"aC",      Color3f(2.9440999183f, 2.2271502925f, 1.9681668794f), Color3f(0.8874329109f, 0.7993216383f, 0.8152862927f)},
            {"Ag",      Color3f(0.1552646489f, 0.1167232965f, 0.1383806959f), Color3f(4.8283433224f, 3.1222459278f, 2.1469504455f)},
            {"Al",      Color3f(1.6574599595f, 0.8803689579f, 0.5212287346f), Color3f(9.2238691996f, 6.2695232477f, 4.8370012281f)},
            {"AlAs",    Color3f(3.6051023902f, 3.2329365777f, 2.2175611545f), Color3f(0.0006670247f, -0.0004999400f, 0.0074261204f)},
            {"AlSb",    Color3f(-0.048522570f, 4.1427547893f, 4.6697691348f), Color3f(-0.0363741915f, 0.0937665154f, 1.3007390124f)},
            {"Au",      Color3f(0.1431189557f, 0.3749570432f, 1.4424785571f), Color3f(3.9831604247f, 2.3857207478f, 1.6032152899f)},
            {"Be",      Color3f(4.1850592788f, 3.1850604423f, 2.7840913457f), Color3f(3.8354398268f, 3.0101260162f, 2.8690088743f)},
            {"Cr",      Color3f(4.3696828663f, 2.9167024892f, 1.6547005413f), Color3f(5.2064337956f, 4.2313645277f, 3.7549467933f)},
            {"CsI",     Color3f(2.1449030413f, 1.7023164587f, 1.6624194173f), Color3f(0.0000000000f, 0.0000000000f, 0.0000000000f)},
            {"Cu",      Color3f(0.2004376970f, 0.9240334304f, 1.1022119527f), Color3f(3.9129485033f, 2.4528477015f, 2.1421879552f)},
            {"Cu2O",    Color3f(3.5492833755f, 2.9520622449f, 2.7369202137f), Color3f(0.1132179294f, 0.1946659670f, 0.6001681264f)},
            {"CuO",     Color3f(3.2453822204f, 2.4496293965f, 2.1974114493f), Color3f(0.5202739621f, 0.5707372756f, 0.7172250613f)},
            {"d-C",     Color3f(2.7112524747f, 2.3185812849f, 2.2288565009f), Color3f(0.0000000000f, 0.0000000000f, 0.0000000000f)},
            {"Hg",      Color3f(2.3989314904f, 1.4400254917f, 0.9095512090f), Color3f(6.3276269444f, 4.3719414152f, 3.4217899270f)},
            {"HgTe",    Color3f(4.7795267752f, 3.2309984581f, 2.6600252401f), Color3f(1.6319827058f, 1.5808189339f, 1.7295753852f)},
            {"Ir",      Color3f(3.0864098394f, 2.0821938440f, 1.6178866805f), Color3f(5.5921510077f, 4.0671757150f, 3.2672611269f)},
            {"K",       Color3f(0.0640493070f, 0.0464100621f, 0.0381842017f), Color3f(2.1042155920f, 1.3489364357f, 0.9132113889f)},
            {"Li",      Color3f(0.2657871942f, 0.1956102432f, 0.2209198538f), Color3f(3.5401743407f, 2.3111306542f, 1.6685930000f)},
            {"MgO",     Color3f(2.0895885542f, 1.6507224525f, 1.5948759692f), Color3f(0.0000000000f, -0.0000000000f, 0.0000000000f)},
            {"Mo",      Color3f(4.4837010280f, 3.5254578255f, 2.7760769438f), Color3f(4.1111307988f, 3.4208716252f, 3.1506031404f)},
            {"Na",      Color3f(0.0602665320f, 0.0561412435f, 0.0619909494f), Color3f(3.1792906496f, 2.1124800781f, 1.5790940266f)},
            {"Nb",      Color3f(3.4201353595f, 2.7901921379f, 2.3955856658f), Color3f(3.4413817900f, 2.7376437930f, 2.5799132708f)},
            {"Ni",      Color3f(2.3672753521f, 1.6633583302f, 1.4670554172f), Color3f(4.4988329911f, 3.0501643957f, 2.3454274399f)},
            {"Rh",      Color3f(2.5857954933f, 1.8601866068f, 1.5544279524f), Color3f(6.7822927110f, 4.7029501026f, 3.9760892461f)},
            {"Se-e",    Color3f(5.7242724833f, 4.1653992967f, 4.0816099264f), Color3f(0.8713747439f, 1.1052845009f, 1.5647788766f)},
            {"Se",      Color3f(4.0592611085f, 2.8426947380f, 2.8207582835f), Color3f(0.7543791750f, 0.6385150558f, 0.5215872029f)},
            {"SiC",     Color3f(3.1723450205f, 2.5259677964f, 2.4793623897f), Color3f(0.0000007284f, -0.0000006859f, 0.0000100150f)},
            {"SnTe",    Color3f(4.5251865890f, 1.9811525984f, 1.2816819226f), Color3f(0.0000000000f, 0.0000000000f, 0.0000000000f)},
            {"Ta",      Color3f(2.0625846607f, 2.3930915569f, 2.6280684948f), Color3f(2.4080467973f, 1.7413705864f, 1.9470377016f)},
            {"Te-e",    Color3f(7.5090397678f, 4.2964603080f, 2.3698732430f), Color3f(5.5842076830f, 4.9476231084f, 3.9975145063f)},
            {"Te",      Color3f(7.3908396088f, 4.4821028985f, 2.6370708478f), Color3f(3.2561412892f, 3.5273908133f, 3.2921683116f)},
            {"ThF4",    Color3f(1.8307187117f, 1.4422274283f, 1.3876488528f), Color3f(0.0000000000f, 0.0000000000f, 0.0000000000f)},
            {"TiC",     Color3f(3.7004673762f, 2.8374356509f, 2.5823030278f), Color3f(3.2656905818f, 2.3515586388f, 2.1727857800f)},
            {"TiN",     Color3f(1.6484691607f, 1.1504482522f, 1.3797795097f), Color3f(3.3684596226f, 1.9434888540f, 1.1020123347f)},
            {"TiO2-e",  Color3f(3.1065574823f, 2.5131551146f, 2.5823844157f), Color3f(0.0000289537f, -0.0000251484f, 0.0001775555f)},
            {"TiO2",    Color3f(3.4566203131f, 2.8017076558f, 2.9051485020f), Color3f(0.0001026662f, -0.0000897534f, 0.0006356902f)},
            {"VC",      Color3f(3.6575665991f, 2.7527298065f, 2.5326814570f), Color3f(3.0683516659f, 2.1986687713f, 1.9631816252f)},
            {"VN",      Color3f(2.8656011588f, 2.1191817791f, 1.9400767149f), Color3f(3.0323264950f, 2.0561075580f, 1.6162930914f)},
            {"V",       Color3f(4.2775126218f, 3.5131538236f, 2.7611257461f), Color3f(3.4911844504f, 2.8893580874f, 3.1116965117f)},
            {"W",       Color3f(4.3707029924f, 3.3002972445f, 2.9982666528f), Color3f(3.5006778591f, 2.6048652781f, 2.2731930614f)},
        };
// the number of materials
static const int materialConst = 40;

/**
 * \brief Extracted the material propertoers
 *
 * \param name
 *     The name of the material
 * \param eta
 *     The refraction coefficient
 * \param k
 *     The absorbtion coefficient
 *
 * \return
 *     If the material was found.
 */
bool getMaterialProperties(const std::string &name, Color3f &eta, Color3f &k){
    for (int i = 0; i < materialConst; ++i) {
        if (matrialList[i].name == name) {
            eta = matrialList[i].eta;
            k   = matrialList[i].k;
            return true;
        }
    }
    return false;
}

NORI_NAMESPACE_END

#endif /* __NORI_METALCONSTS_H */

