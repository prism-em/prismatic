#include "prism_colormapper.h"
#include <QColor>
#include <cmath>
#include <iostream>
#include <map>
namespace Prismatic
{

void Colormapper::setColormap(const Colormap &cmap) { this->colormap = cmap; }

QRgb Colormapper::getColor(const double value, const double contrastMin, const double contrastMax)
{

    Color<double> color_fraction;
    Color<unsigned char> c;
    //        std::cout << "value  " << value << std::endl;
    //        std::cout << "contrastMax  " << contrastMax << std::endl;
    //        std::cout << "contrastMin  " << contrastMin << std::endl;
    if (value <= contrastMin)
    {
        // use first color in colormap
        color_fraction = this->colormap[0];

        // convert to uint8
        c = Color<unsigned char>{(unsigned char)(color_fraction.r * 255.0),
                                 (unsigned char)(color_fraction.g * 255.0),
                                 (unsigned char)(color_fraction.b * 255.0)};
    }
    else if (value >= contrastMax)
    {
        // use last color in colormap
        color_fraction = this->colormap[this->colormap.size() - 1];

        // convert to uint8
        c = Color<unsigned char>{(unsigned char)(color_fraction.r * 255.0),
                                 (unsigned char)(color_fraction.g * 255.0),
                                 (unsigned char)(color_fraction.b * 255.0)};
    }
    else
    {
        // linearly interpolate the color value based upon the two nearest positional neighbors

        // determine what fraction of max value the value is
        float fractional_pos = (value - contrastMin) / (contrastMax - contrastMin);

        // map this fraction to a color (or decimal index between two colors)
        float colorIndex = fractional_pos * (this->colormap.size() - 1);

        // determine the two relevant colors
        int lowerColorIndex = (int)std::floor(colorIndex);
        int upperColorIndex = (int)std::ceil(colorIndex);
        Color<double> cLow, cHigh;
        cLow = this->colormap[lowerColorIndex];
        cHigh = this->colormap[upperColorIndex];

        // determine the weights
        double weight1, weight2;
        weight2 = colorIndex - (float)lowerColorIndex;
        weight1 = 1 - weight2;
        //            std::cout << "weight1 = " << weight1 << std::endl;
        //            std::cout << "weight2 = " << weight2 << std::endl;
        //            std::cout << "lowerColorIndex = " << lowerColorIndex << std::endl;
        //            std::cout << "upperColorIndex = " << upperColorIndex << std::endl;

        // make the interpolated color
        c = Color<unsigned char>{(unsigned char)((weight1 * cLow.r + weight2 * cHigh.r) * 255.0),
                                 (unsigned char)((weight1 * cLow.g + weight2 * cHigh.g) * 255.0),
                                 (unsigned char)((weight1 * cLow.b + weight2 * cHigh.b) * 255.0)};
    }
    return qRgba(c.r, c.g, c.b, 255);
    //          return qRgba(0, 0, 0, 255);
}

} // namespace Prismatic
