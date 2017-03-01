//
// Created by AJ Pryor on 2/21/17.
//
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <map>
#include <random>
#include <thread>
#include "emdSTEM.h"
#include "ArrayND.h"
#include "projPot.h"
using namespace std;
namespace PRISM {
	template <class T>
	using Array3D = ArrayND<3, std::vector<T> >;
	template <class T>
	using Array2D = ArrayND<2, std::vector<T> >;
	template <class T>
	using Array1D = ArrayND<1, std::vector<T> >;

	template <class T>
	void fetch_potentials(PRISM::ArrayND<3, std::vector<T> >& potentials,
	                      const vector<size_t>& atomic_species,
	                      const ArrayND<1, std::vector<T> >& xr,
	                      const ArrayND<1, std::vector<T> >& yr){
		ArrayND<2, std::vector<T> > cur_pot;
		for (auto k =0; k < potentials.get_dimk(); ++k){
			cout << "species = " << atomic_species[k] << endl;
			ArrayND<2, std::vector<T> > cur_pot = projPot(atomic_species[k], xr, yr);
			for (auto j = 0; j < potentials.get_dimj(); ++j){
				for (auto i = 0; i < potentials.get_dimi(); ++i){
					potentials.at(k,j,i) = cur_pot.at(j,i);
				}
			}
		}
	}

	template <class T>
	vector<size_t> get_unique_atomic_species(emdSTEM<T>& pars){
		// helper function to get the unique atomic species
		vector<size_t> unique_atoms = vector<size_t>(pars.atoms.size(),0);
		for (auto i = 0; i < pars.atoms.size(); ++i)unique_atoms[i] = pars.atoms[i].species;
		sort(unique_atoms.begin(), unique_atoms.end());
		vector<size_t>::iterator it = unique(unique_atoms.begin(), unique_atoms.end());
		unique_atoms.resize(distance(unique_atoms.begin(),it));
		return unique_atoms;
	}

	template <class T>
	void generateProjectedPotentials(emdSTEM<T>& pars,
                                           const Array3D<T>& potLookup,
                                           const vector<size_t>& unique_species,
                                           const ArrayND<1, vector<long> >& xvec,
                                           const ArrayND<1, vector<long> >& yvec,
                                           const Array1D<T>& uLookup ){
		// splits the atomic coordinates into slices and computes the projected potential for each.
		Array1D<T> x  = zeros_ND<1, T>({pars.atoms.size()});
		Array1D<T> y  = zeros_ND<1, T>({pars.atoms.size()});
		Array1D<T> z  = zeros_ND<1, T>({pars.atoms.size()});
		Array1D<T> ID = zeros_ND<1, T>({pars.atoms.size()});

		for (auto i = 0; i < pars.atoms.size(); ++i){
			x[i]  = pars.atoms[i].x * pars.cellDim[2];
			y[i]  = pars.atoms[i].y * pars.cellDim[1];
			z[i]  = pars.atoms[i].z * pars.cellDim[0];
			ID[i] = pars.atoms[i].species;
		}

		// compute the z-slice index for each atom
		auto max_z = std::max_element(z.begin(), z.end());
		Array1D<T> zPlane(z);
		std::transform(zPlane.begin(), zPlane.end(), zPlane.begin(), [&max_z, &pars](T &t_z) {
			return round((-t_z + *max_z) / pars.sliceThickness + 0.5) - 1; // If the +0.5 was to make the first slice z=1 not 0, can drop the +0.5 and -1
		});
		max_z = std::max_element(zPlane.begin(), zPlane.end());
		pars.numPlanes = *max_z + 1;

		pars.pot = zeros_ND<3, T>({pars.numPlanes,pars.imageSize[0], pars.imageSize[1] });


		// create a key-value map to match the atomic Z numbers with their place in the potential lookup table
		map<size_t, size_t> Z_lookup;
		for (auto i = 0; i < unique_species.size(); ++i)Z_lookup[unique_species[i]] = i;
		cout << "Z_lookup[0]  = " << Z_lookup[0] << endl;
		cout << "Z_lookup[78]  = " << Z_lookup[78] << endl;
		cout << "Z_lookup[6]  = " << Z_lookup[6] <<endl;


		// loop over each plane, perturb the atomic positions, and place the corresponding potential at each location

		// parallel calculation of each individual slice
		std::vector<std::thread> workers;
		workers.reserve(pars.numPlanes);
		for (long a0 = 0; a0 < pars.numPlanes; ++a0){
			workers.emplace_back(thread([&pars, &x, &y, &z, &ID, &Z_lookup, &xvec, &zPlane, &yvec,&potLookup,&uLookup,a0](){
				std::default_random_engine de(time(0));
				normal_distribution<T> randn(0,1);
				ArrayND<1, vector<long> > xp;
				ArrayND<1, vector<long> > yp;
				Array2D<T> projPot = zeros_ND<2, T>({pars.imageSize[0], pars.imageSize[1]});
				for (auto a2 = 0; a2 < x.size(); ++a2){
					if (zPlane[a2]==a0){
						const size_t cur_Z = Z_lookup[ID[a2]];
						const T X = round((x[a2]) / pars.pixelSize[1]);
//						const T X = round((x[a2] + randn(de)*uLookup[cur_Z]) / pars.pixelSize[0]);
						xp = xvec + (long)X;
						for (auto& i:xp)i%=pars.imageSize[1];
						const T Y = round((y[a2])/pars.pixelSize[0]);
//						const T Y = round((y[a2] + randn(de)*uLookup[cur_Z]) / pars.pixelSize[1]);
						yp = yvec + (long)Y;
						for (auto& i:yp)i%=pars.imageSize[0];
						for (auto ii = 0; ii < xp.size(); ++ii){
							for (auto jj = 0; jj < yp.size(); ++jj){
								projPot.at(yp[jj],xp[ii]) += potLookup.at(cur_Z,jj,ii);
							}
						}
					}
				}
				copy(projPot.begin(), projPot.end(),&pars.pot.at(a0,0,0));

			}));

		}
		for (auto &t:workers)t.join();


		//serial version
//				Array2D<T> projPot = zeros_ND<2, T>({pars.imageSize[0], pars.imageSize[1]});
//		std::default_random_engine de(time(0));
//				normal_distribution<T> randn(0,1);
//				ArrayND<1, vector<long> > xp;
//				ArrayND<1, vector<long> > yp;
//		for (auto a0 = 0; a0 < pars.numPlanes; ++a0){
//			memset((void*)&projPot[0], 0, projPot.size() * sizeof(T));
//			for (auto a2 = 0; a2 < x.size(); ++a2){
//				if (zPlane[a2]==a0){
//					const size_t cur_Z = Z_lookup[ID[a2]];
//					const T X = round((x[a2] + randn(de)*uLookup[cur_Z]) / pars.pixelSize[0]);
//					xp = xvec + (long)X;
//					for (auto& i:xp)i%=pars.imageSize[0];
//					const T Y = round((y[a2] + randn(de)*uLookup[cur_Z]) / pars.pixelSize[1]);
//					yp = yvec + (long)Y;
//					for (auto& i:yp)i%=pars.imageSize[1];
//					for (auto ii = 0; ii < xp.size(); ++ii){
//						for (auto jj = 0; jj < yp.size(); ++jj){
//							projPot.at(yp[jj],xp[ii]) += potLookup.at(cur_Z,jj,ii);
//						}
//					}
//					//const T X = x[a2] + randn(de)*uLookup[Z_lookup[ID[a2]]];
////					transform(xvec.begin(), xvec.end(), xvec.begin(),[](){
////
////						;});
//				}
//			}
//			copy(projPot.begin(), projPot.end(),&pars.pot.at(a0,0,0));
//		}

		cout << "pars.numPlans = " << pars.numPlanes << endl;
		cout << "zPlane[1] = " << zPlane[1] << endl;
		cout << "zPlane[2] = " << zPlane[2] << endl;
		cout << "zPlane[3] = " << zPlane[3] << endl;
		cout << "zPlane[100] = " << zPlane[100] << endl;
		cout << "zPlane[1000] = " << zPlane[1000] << endl;

		cout << "x[1] = " << x[1] << endl;
		cout << "y[1] = " << y[1] << endl;
		cout << "z[1] = " << z[1] << endl;
//		return projPot;
	};

	template <class T>
	void PRISM01(emdSTEM<T>& pars){
		//builds atomic potentials

		using Array3D = ArrayND<3, std::vector<T> >;
		using Array2D = ArrayND<2, std::vector<T> >;
		using Array1D = ArrayND<1, std::vector<T> >;

		//forward declare helper functions

//		void fetch_potentials(PRISM::ArrayND<3, std::vector<T> >& potentials,
//                                     const vector<size_t>& atomic_species,
//                                     const ArrayND<1, std::vector<T> >& xr,
//                                     const ArrayND<1, std::vector<T> >& yr);
		ArrayND<2, std::vector<T> > projPot(const size_t&, const Array1D&, const Array1D&);

		cout << "Entering PRISM01" << endl;
		T yleng = std::ceil(pars.potBound / pars.pixelSize[0]);
		T xleng = std::ceil(pars.potBound / pars.pixelSize[1]);
		ArrayND<1, vector<long> > xvec(vector<long>(2*(size_t)xleng + 1, 0),{2*(size_t)xleng + 1});
		ArrayND<1, vector<long> > yvec(vector<long>(2*(size_t)yleng + 1, 0),{2*(size_t)yleng + 1});
		{
			T tmpx = -xleng;
			T tmpy = -yleng;
			for (auto &i : xvec)i = tmpx++;
			for (auto &j : yvec)j = tmpy++;
		}
		Array1D xr(vector<T>(2*(size_t)xleng + 1, 0),{2*(size_t)xleng + 1});
		Array1D yr(vector<T>(2*(size_t)yleng + 1, 0),{2*(size_t)yleng + 1});
		for (auto i=0; i < xr.size(); ++i)xr[i] = (T)xvec[i] * pars.pixelSize[1];
		for (auto j=0; j < yr.size(); ++j)yr[j] = (T)yvec[j] * pars.pixelSize[0];

		vector<size_t> unique_species = get_unique_atomic_species(pars);
		cout << "unique_species[0] = " << unique_species[0] << endl;
		cout << "unique_species[1] = " << unique_species[1] << endl;

		Array1D uLookup   = zeros_ND<1, T>({unique_species.size()});
		for (auto i = 0; i < unique_species.size(); ++i){
			uLookup[i] = pars.u[unique_species[i]] > 0 ? pars.u[unique_species[i]] : 0.05;
		}
		Array3D potLookup = zeros_ND<3, T>({unique_species.size(), 2*(size_t)yleng + 1, 2*(size_t)xleng + 1});
		fetch_potentials(potLookup, unique_species, xr, yr);

		generateProjectedPotentials(pars, potLookup, unique_species, xvec, yvec, uLookup);

#ifndef NDEBUG
//		for (auto& i :uLookup)cout<<i << endl;
#endif //NDEBUG

	}
}
