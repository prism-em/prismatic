//
// Created by AJ Pryor on 2/21/17.
//
#include <iostream>
#include <algorithm>
#include "emdSTEM.h"
#include "ArrayND.h"
using namespace std;
namespace PRISM {
	template <class T>
	vector<size_t> get_unique_atomic_species(emdSTEM<T>&);

	template <class T>
	void PRISM01(emdSTEM<T>& pars){
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
		ArrayND<1, vector<T> > xr(vector<T>(2*(size_t)xleng + 1, 0),{2*(size_t)xleng + 1});
		ArrayND<1, vector<T> > yr(vector<T>(2*(size_t)yleng + 1, 0),{2*(size_t)yleng + 1});
		for (auto i=0; i < xr.size(); ++i)xr[i] = (T)xvec[i] * pars.pixelSize[0];
		for (auto j=0; j < yr.size(); ++j)yr[j] = (T)yvec[j] * pars.pixelSize[1];

		vector<size_t> unique_species = get_unique_atomic_species(pars);

#ifndef NDEBUG
		cout << "pars.pixelSize[0] = " << pars.pixelSize[0] << endl;
		cout << "pars.pixelSize[0] = " << pars.pixelSize[1] << endl;
		//for (auto& j : yvec)cout<<j<<endl;
		for (auto& j : yr)cout<<j<<endl;
		cout<<yr.size()<<endl;
		cout<<xr.size()<<endl;
		cout<<yvec.size()<<endl;
		for (auto& i:unique_species)cout<<i<<endl;
#endif //NDEBUG

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
}
