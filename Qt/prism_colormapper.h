#ifndef PRISMCOLORMAPS_H
#define PRISMCOLORMAPS_H
#include <vector>


namespace Prismatic{
	struct Color{unsigned char r,g,b;};
	using Colormap = std::vector<Color>;

	const Colormap nullColormap = {{0, 0, 0}};
	const Colormap GrayscaleColormap = {{0,0,0}, {1,1,1}};

	class Colormapper{
	public:
		Colormapper(const Colormap& cmap) : colormap(cmap){};
	    Color getColor(const float value, const float contrastMax, const float contrastMin);
	    void setColormap(const Colormap& cmap);
	private:
		Colormap colormap;
	};
}

#endif //PRISMCOLORMAPS_H