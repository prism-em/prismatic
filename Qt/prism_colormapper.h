#ifndef PRISMCOLORMAPS_H
#define PRISMCOLORMAPS_H
#include <vector>
#include <QColor>

namespace Prismatic{

    template<class T>
    struct Color{T r,g,b;};
    using Colormap = std::vector< Color<double> >;

    const Colormap nullColormap      = {{0.0, 0.0, 0.0}};
    const Colormap GrayscaleColormap = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    const Colormap RedColormap       = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    const Colormap BlueColormap      = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    const Colormap GreenColormap     = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    const Colormap MixColormap     = {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}};

	class Colormapper{
	public:
		Colormapper(const Colormap& cmap) : colormap(cmap){};
        QRgb getColor(const double value, const double contrastMin, const double contrastMax);
	    void setColormap(const Colormap& cmap);
	private:
		Colormap colormap;
	};
}

#endif //PRISMCOLORMAPS_H
