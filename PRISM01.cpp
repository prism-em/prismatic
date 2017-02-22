//
// Created by AJ Pryor on 2/21/17.
//
#include <iostream>
#include <algorithm>
#include "emdSTEM.h"
#include "ArrayND.h"
#include "projPot.h"
using namespace std;
namespace PRISM {
	template <class T>
	void fetch_potentials(PRISM::ArrayND<3, std::vector<T> >& potentials,
	                      const vector<size_t>& atomic_species,
	                      const ArrayND<1, std::vector<T> >& xr,
	                      const ArrayND<1, std::vector<T> >& yr){
		ArrayND<2, std::vector<T> > cur_pot;
		for (auto a0 =0; a0 < potentials.get_nrows(); ++a0){
			ArrayND<2, std::vector<T> > cur_pot = projPot(atomic_species[a0], xr, yr);
			for (auto j = 0; j < potentials.get_ncols(); ++j){
				for (auto k = 0; k < potentials.get_nlayers(); ++k){
					potentials.at(a0,j,k) = (T)(rand() % 100) / 10;
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
		T xleng = std::ceil(pars.potBound / pars.pixelSize[0]);
		T yleng = std::ceil(pars.potBound / pars.pixelSize[1]);
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
		for (auto i=0; i < xr.size(); ++i)xr[i] = (T)xvec[i] * pars.pixelSize[0];
		for (auto j=0; j < yr.size(); ++j)yr[j] = (T)yvec[j] * pars.pixelSize[1];

		vector<size_t> unique_species = get_unique_atomic_species(pars);


		Array1D uLookup   = zeros_ND<1, T>({unique_species.size()});
		for (auto i = 0; i < unique_species.size(); ++i){
			uLookup[i] = pars.u[unique_species[i]] > 0 ? pars.u[unique_species[i]] : 0.05;
		}
		Array3D potLookup = zeros_ND<3, T>({unique_species.size(), 2*(size_t)yleng + 1, 2*(size_t)xleng + 1});
		fetch_potentials(potLookup, unique_species, xr, yr);

#ifndef NDEBUG
//		cout << "pars.pixelSize[0] = " << pars.pixelSize[0] << endl;
//		cout << "pars.pixelSize[0] = " << pars.pixelSize[1] << endl;
//		//for (auto& j : yvec)cout<<j<<endl;
//		for (auto& j : yr)cout<<j<<endl;
//		cout<<yr.size()<<endl;
//		cout<<xr.size()<<endl;
//		cout<<yvec.size()<<endl;
//		for (auto& i:unique_species)cout<<i<<endl;
//		cout <<"number of unique atomic species = " << unique_species.size() << endl;
//		cout <<"uLookup.size() = " << uLookup.size() << endl;
//		for (auto& i :potLookup)cout<<i << endl;
//		for (auto& i :uLookup)cout<<i << endl;
#endif //NDEBUG

	}




}
