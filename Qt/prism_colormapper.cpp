#include "prism_colormapper.h"

namespace Prismatic{
	
	void Colormapper::setColormap(const Colormap& cmap){this->colormap = cmap;}

	Color Colormapper::getColor(const float value, const float contrastMax, const float contrastMin) {
		
		if (value < contrastMin) return {0, 0, 0};
		if (value > contrastMax) return {255, 255, 255};
		unsigned char shade =  (unsigned char)( (value - contrastMin) / (contrastMax - contrastMin) * 255);
		return {shade, shade, shade};
	}
	
}
