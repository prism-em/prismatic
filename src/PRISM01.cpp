// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "PRISM01.h"
#include "params.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <map>
#include <vector>
#include <random>
#include <thread>
#include "params.h"
#include "ArrayND.h"
#include "projectedPotential.h"
#include "WorkDispatcher.h"

#ifdef PRISM_BUILDING_GUI
#include "prism_progressbar.h"
#endif


namespace PRISM {
	using namespace std;
	void fetch_potentials(Array3D<PRISM_FLOAT_PRECISION>& potentials,
	                      const vector<size_t>& atomic_species,
	                      const Array1D<PRISM_FLOAT_PRECISION>& xr,
	                      const Array1D<PRISM_FLOAT_PRECISION>& yr){
		Array2D<PRISM_FLOAT_PRECISION> cur_pot;
		for (auto k =0; k < potentials.get_dimk(); ++k){
			cout <<"atomic_species[" <<k<<"] = " << atomic_species[k] << endl;
			Array2D<PRISM_FLOAT_PRECISION> cur_pot = projPot(atomic_species[k], xr, yr);
			for (auto j = 0; j < potentials.get_dimj(); ++j){
				for (auto i = 0; i < potentials.get_dimi(); ++i){
					potentials.at(k,j,i) = cur_pot.at(j,i);
				}
			}
		}
	}

	vector<size_t> get_unique_atomic_species(Parameters<PRISM_FLOAT_PRECISION>& pars){
		// helper function to get the unique atomic species
		vector<size_t> unique_atoms = vector<size_t>(pars.atoms.size(),0);
		for (auto i = 0; i < pars.atoms.size(); ++i)unique_atoms[i] = pars.atoms[i].species;
		sort(unique_atoms.begin(), unique_atoms.end());
		vector<size_t>::iterator it = unique(unique_atoms.begin(), unique_atoms.end());
		unique_atoms.resize(distance(unique_atoms.begin(),it));
		return unique_atoms;
	}

	void generateProjectedPotentials(Parameters<PRISM_FLOAT_PRECISION>& pars,
	                                 const Array3D<PRISM_FLOAT_PRECISION>& potentialLookup,
	                                 const vector<size_t>& unique_species,
	                                 const Array1D<long>& xvec,
	                                 const Array1D<long>& yvec,
	                                 const Array1D<PRISM_FLOAT_PRECISION>& uLookup ){
		// splits the atomic coordinates into slices and computes the projected potential for each.

		// create arrays for the coordinates
		Array1D<PRISM_FLOAT_PRECISION> x  = zeros_ND<1, PRISM_FLOAT_PRECISION>({{pars.atoms.size()}});
		Array1D<PRISM_FLOAT_PRECISION> y  = zeros_ND<1, PRISM_FLOAT_PRECISION>({{pars.atoms.size()}});
		Array1D<PRISM_FLOAT_PRECISION> z  = zeros_ND<1, PRISM_FLOAT_PRECISION>({{pars.atoms.size()}});
		Array1D<PRISM_FLOAT_PRECISION> ID = zeros_ND<1, PRISM_FLOAT_PRECISION>({{pars.atoms.size()}});

		// populate arrays from the atoms structure
		for (auto i = 0; i < pars.atoms.size(); ++i){
			x[i]  = pars.atoms[i].x * pars.meta.cellDim[2];
			y[i]  = pars.atoms[i].y * pars.meta.cellDim[1];
			z[i]  = pars.atoms[i].z * pars.meta.cellDim[0];
			ID[i] = pars.atoms[i].species;
		}


		// compute the z-slice index for each atom
		auto max_z = std::max_element(z.begin(), z.end());
		Array1D<PRISM_FLOAT_PRECISION> zPlane(z);
		std::transform(zPlane.begin(), zPlane.end(), zPlane.begin(), [&max_z, &pars](PRISM_FLOAT_PRECISION &t_z) {
			return round((-t_z + *max_z) / pars.meta.sliceThickness + 0.5) - 1; // If the +0.5 was to make the first slice z=1 not 0, can drop the +0.5 and -1
		});
		max_z = std::max_element(zPlane.begin(), zPlane.end());
		pars.numPlanes = *max_z + 1;
#ifdef PRISM_BUILDING_GUI
		pars.progressbar->signalPotentialUpdate(0, pars.numPlanes);
#endif

		// initialize the potential array
		pars.pot = zeros_ND<3, PRISM_FLOAT_PRECISION>({{pars.numPlanes, pars.imageSize[0], pars.imageSize[1]}});

		// create a key-value map to match the atomic Z numbers with their place in the potential lookup table
		map<size_t, size_t> Z_lookup;
		for (auto i = 0; i < unique_species.size(); ++i)Z_lookup[unique_species[i]] = i;

		//loop over each plane, perturb the atomic positions, and place the corresponding potential at each location
		// using parallel calculation of each individual slice
		std::vector<std::thread> workers;
		workers.reserve(pars.meta.NUM_THREADS);

//		setWorkStartStop(0, pars.numPlanes, 1);
		WorkDispatcher dispatcher(0, pars.numPlanes, 1);
		for (long t = 0; t < pars.meta.NUM_THREADS; ++t){
			cout << "Launching thread #" << t << " to compute projected potential slices\n";
			workers.push_back(thread([&pars, &x, &y, &z, &ID, &Z_lookup, &xvec,
											 &zPlane, &yvec,&potentialLookup,&uLookup, &dispatcher](){
				// create a random number generator to simulate thermal effects
//				std::default_random_engine de(time(0));
				cout <<"pars.meta.random_seed = " << pars.meta.random_seed<< endl;
				std::default_random_engine de(pars.meta.random_seed);
				normal_distribution<PRISM_FLOAT_PRECISION> randn(0,1);
				Array1D<long> xp;
				Array1D<long> yp;

				size_t currentBeam, stop;
                currentBeam=stop=0;
				//while (getWorkID(pars, currentBeam, stop)) { // synchronously get work assignment
				while (dispatcher.getWork(currentBeam, stop)) { // synchronously get work assignment
					Array2D<PRISM_FLOAT_PRECISION> projectedPotential = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.imageSize[0], pars.imageSize[1]}});
					while (currentBeam != stop) {
						for (auto a2 = 0; a2 < x.size(); ++a2) {
							if (zPlane[a2] == currentBeam) {
								const long dim0 = (long) pars.imageSize[0];
								const long dim1 = (long) pars.imageSize[1];
								const size_t cur_Z = Z_lookup[ID[a2]];
								PRISM_FLOAT_PRECISION X, Y;
								if (pars.meta.include_thermal_effects) {
									X = round(
											(x[a2] + randn(de) * uLookup[cur_Z]) / pars.pixelSize[0]);
									Y = round(
											(y[a2] + randn(de) * uLookup[cur_Z]) / pars.pixelSize[1]);
								} else {
									X = round((x[a2]) / pars.pixelSize[1]); // this line uses no thermal factor
									Y = round((y[a2]) / pars.pixelSize[0]); // this line uses no thermal factor
								}
								xp = xvec + (long) X;
								for (auto &i:xp)i = (i % dim1 + dim1) % dim1; // make sure to get a positive value

								yp = yvec + (long) Y;
								for (auto &i:yp) i = (i % dim0 + dim0) % dim0;// make sure to get a positive value
								for (auto ii = 0; ii < xp.size(); ++ii) {
									for (auto jj = 0; jj < yp.size(); ++jj) {
										projectedPotential.at(yp[jj], xp[ii]) += potentialLookup.at(cur_Z, jj, ii);
									}
								}
							}
						}
						copy(projectedPotential.begin(), projectedPotential.end(),&pars.pot.at(currentBeam,0,0));
#ifdef PRISM_BUILDING_GUI
                        pars.progressbar->signalPotentialUpdate(currentBeam, pars.numPlanes);
#endif //PRISM_BUILDING_GUI
						++currentBeam;
					}
				}

			}));

		}
		cout << "Waiting for threads...\n";
		for (auto &t:workers)t.join();
#ifdef PRISM_BUILDING_GUI
		pars.progressbar->setProgress(100);
#endif //PRISM_BUILDING_GUI
	};

	void PRISM01(Parameters<PRISM_FLOAT_PRECISION>& pars){

		//builds atomic potentials

		// setup some coordinates
		cout << "Entering PRISM01" << endl;
		PRISM_FLOAT_PRECISION yleng = std::ceil(pars.meta.potBound / pars.pixelSize[0]);
		PRISM_FLOAT_PRECISION xleng = std::ceil(pars.meta.potBound / pars.pixelSize[1]);
		ArrayND<1, vector<long> > xvec(vector<long>(2*(size_t)xleng + 1, 0),{{2*(size_t)xleng + 1}});
		ArrayND<1, vector<long> > yvec(vector<long>(2*(size_t)yleng + 1, 0),{{2*(size_t)yleng + 1}});
		{
			PRISM_FLOAT_PRECISION tmpx = -xleng;
			PRISM_FLOAT_PRECISION tmpy = -yleng;
			for (auto &i : xvec)i = tmpx++;
			for (auto &j : yvec)j = tmpy++;
		}
		Array1D<PRISM_FLOAT_PRECISION> xr(vector<PRISM_FLOAT_PRECISION>(2*(size_t)xleng + 1, 0),{{2*(size_t)xleng + 1}});
		Array1D<PRISM_FLOAT_PRECISION> yr(vector<PRISM_FLOAT_PRECISION>(2*(size_t)yleng + 1, 0),{{2*(size_t)yleng + 1}});
		for (auto i=0; i < xr.size(); ++i)xr[i] = (PRISM_FLOAT_PRECISION)xvec[i] * pars.pixelSize[1];
		for (auto j=0; j < yr.size(); ++j)yr[j] = (PRISM_FLOAT_PRECISION)yvec[j] * pars.pixelSize[0];

		vector<size_t> unique_species = get_unique_atomic_species(pars);
		Array1D<PRISM_FLOAT_PRECISION> uLookup   = zeros_ND<1, PRISM_FLOAT_PRECISION>({{unique_species.size()}});
		for (auto i = 0; i < unique_species.size(); ++i){
			uLookup[i] = pars.u[unique_species[i]] > 0 ? pars.u[unique_species[i]] : 0.05;
		}
		
		cout <<"lookups\n";
		// initialize the lookup table
		Array3D<PRISM_FLOAT_PRECISION> potentialLookup = zeros_ND<3, PRISM_FLOAT_PRECISION>({{unique_species.size(), 2*(size_t)yleng + 1, 2*(size_t)xleng + 1}});
		
		cout <<"fetch\n";
		// precompute the unique potentials
		fetch_potentials(potentialLookup, unique_species, xr, yr);
		cout <<"project\n";

		// populate the slices with the projected potentials
		generateProjectedPotentials(pars, potentialLookup, unique_species, xvec, yvec, uLookup);
	}
}
