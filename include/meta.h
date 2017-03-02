//
// Created by AJ Pryor on 3/2/17.
//

#ifndef PRISM_META_H
#define PRISM_META_H
namespace PRISM{
	template <class T>
	class Metadata{
	public:
		Metadata(){
			interpolationFactor = 5;
			filename_atoms      = "";
			filename_output     = "";
		}
		size_t interpolationFactor;
		std::string filename_atoms;
		std::string filename_output;
	};
}
#endif //PRISM_META_H
