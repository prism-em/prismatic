// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#include "PRISM01_calcPotential.h"
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
#include "utility.h"
#include "fileIO.h"
#include "fftw3.h"
#include <complex>

#ifdef PRISMATIC_BUILDING_GUI
#include "prism_progressbar.h"
#endif

namespace Prismatic
{

using namespace std;
mutex potentialWriteLock;
extern mutex fftw_plan_lock;

void fetch_potentials(Array3D<PRISMATIC_FLOAT_PRECISION> &potentials,
					  const vector<size_t> &atomic_species,
					  const Array1D<PRISMATIC_FLOAT_PRECISION> &xr,
					  const Array1D<PRISMATIC_FLOAT_PRECISION> &yr)
{
	Array2D<PRISMATIC_FLOAT_PRECISION> cur_pot;
	for (auto k = 0; k < potentials.get_dimk(); ++k)
	{
		Array2D<PRISMATIC_FLOAT_PRECISION> cur_pot = projPot(atomic_species[k], xr, yr);
		for (auto j = 0; j < potentials.get_dimj(); ++j)
		{
			for (auto i = 0; i < potentials.get_dimi(); ++i)
			{
				potentials.at(k, j, i) = cur_pot.at(j, i);
			}
		}
	}
}

void fetch_potentials3D(Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> &potentials,
					  const vector<size_t> &atomic_species,
					  const Array1D<PRISMATIC_FLOAT_PRECISION> &xr,
					  const Array1D<PRISMATIC_FLOAT_PRECISION> &yr,
					  const Array1D<PRISMATIC_FLOAT_PRECISION> &zr)
{
	Array3D<PRISMATIC_FLOAT_PRECISION> cur_pot;
	PRISMATIC_FFTW_INIT_THREADS();
	for (auto l = 0; l < potentials.get_diml(); l++)
	{
		Array3D<PRISMATIC_FLOAT_PRECISION> cur_pot = kirklandPotential3D(atomic_species[l], xr, yr, zr);
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> fstore = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{cur_pot.get_dimj(), cur_pot.get_dimi()}});
		for (auto k = 0; k < cur_pot.get_dimk(); k++)
		{
			//fourier transform potentials in K loop since we only transform in x, y
			for(auto j = 0; j < cur_pot.get_dimj(); j ++)
			{
				for(auto i = 0; i < cur_pot.get_dimi(); i++)
				{
					fstore.at(j,i).real(cur_pot.at(k,j,i));
				}
			}
			unique_lock<mutex> gatekeeper(fftw_plan_lock);
			PRISMATIC_FFTW_PLAN plan_forward = PRISMATIC_FFTW_PLAN_DFT_2D(cur_pot.get_dimj(), cur_pot.get_dimi(),
																	reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&fstore[0]),
																	reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&potentials.at(l,k,0,0)),
																	FFTW_FORWARD,
																	FFTW_ESTIMATE);

			gatekeeper.unlock();
			PRISMATIC_FFTW_EXECUTE(plan_forward);

			gatekeeper.lock();
			PRISMATIC_FFTW_DESTROY_PLAN(plan_forward);
		}
	}
	PRISMATIC_FFTW_CLEANUP_THREADS();
}

vector<size_t> get_unique_atomic_species(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// helper function to get the unique atomic species
	vector<size_t> unique_atoms = vector<size_t>(pars.atoms.size(), 0);
	for (auto i = 0; i < pars.atoms.size(); ++i)
		unique_atoms[i] = pars.atoms[i].species;
	sort(unique_atoms.begin(), unique_atoms.end());
	vector<size_t>::iterator it = unique(unique_atoms.begin(), unique_atoms.end());
	unique_atoms.resize(distance(unique_atoms.begin(), it));
	return unique_atoms;
}

void generateProjectedPotentials(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
								 const Array3D<PRISMATIC_FLOAT_PRECISION> &potentialLookup,
								 const vector<size_t> &unique_species,
								 const Array1D<long> &xvec,
								 const Array1D<long> &yvec)
{
	// splits the atomic coordinates into slices and computes the projected potential for each.

	// create arrays for the coordinates
	Array1D<PRISMATIC_FLOAT_PRECISION> x = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> y = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> z = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> ID = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> sigma = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> occ = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});

	// populate arrays from the atoms structure
	for (auto i = 0; i < pars.atoms.size(); ++i)
	{
		x[i] = pars.atoms[i].x * pars.tiledCellDim[2];
		y[i] = pars.atoms[i].y * pars.tiledCellDim[1];
		z[i] = pars.atoms[i].z * pars.tiledCellDim[0];
		ID[i] = pars.atoms[i].species;
		sigma[i] = pars.atoms[i].sigma;
		occ[i] = pars.atoms[i].occ;
	}

	// compute the z-slice index for each atom
	Array1D<PRISMATIC_FLOAT_PRECISION> zPlane(z);
	std::transform(zPlane.begin(), zPlane.end(), zPlane.begin(), [&pars](PRISMATIC_FLOAT_PRECISION &t_z) {
		return round((-t_z + pars.tiledCellDim[0]) / pars.meta.sliceThickness + 0.5) - 1; // If the +0.5 was to make the first slice z=1 not 0, can drop the +0.5 and -1
	});
	auto max_z = std::max_element(zPlane.begin(), zPlane.end());
	pars.numPlanes = *max_z + 1;

	//check if intermediate output was specified, if so, create index of output slices
	if (pars.meta.numSlices == 0)
	{
		pars.numSlices = pars.numPlanes;
	}

#ifdef PRISMATIC_BUILDING_GUI
	pars.progressbar->signalPotentialUpdate(0, pars.numPlanes);
#endif

	// initialize the potential array
	pars.pot = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{pars.numPlanes, pars.imageSize[0], pars.imageSize[1]}});

	// create a key-value map to match the atomic Z numbers with their place in the potential lookup table
	map<size_t, size_t> Z_lookup;
	for (auto i = 0; i < unique_species.size(); ++i)
		Z_lookup[unique_species[i]] = i;

	//loop over each plane, perturb the atomic positions, and place the corresponding potential at each location
	// using parallel calculation of each individual slice
	std::vector<std::thread> workers;
	workers.reserve(pars.meta.numThreads);

	WorkDispatcher dispatcher(0, pars.numPlanes);
	for (long t = 0; t < pars.meta.numThreads; ++t)
	{
		cout << "Launching thread #" << t << " to compute projected potential slices\n";
		workers.push_back(thread([&pars, &x, &y, &z, &ID, &Z_lookup, &xvec, &sigma, &occ,
								  &zPlane, &yvec, &potentialLookup, &dispatcher]()
		{
			// create a random number generator to simulate thermal effects
			// std::cout<<"random seed = " << pars.meta.randomSeed << std::endl;
			// srand(pars.meta.randomSeed);
			// std::default_random_engine de(pars.meta.randomSeed);
			// normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0,1);
			Array1D<long> xp;
			Array1D<long> yp;

			size_t currentSlice, stop;
			currentSlice = stop = 0;
			while (dispatcher.getWork(currentSlice, stop))
			{ // synchronously get work assignment
				Array2D<PRISMATIC_FLOAT_PRECISION> projectedPotential = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[0], pars.imageSize[1]}});
				const long dim0 = (long)pars.imageSize[0];
				const long dim1 = (long)pars.imageSize[1];
				while (currentSlice != stop)
				{

					// create a random number generator to simulate thermal effects
					std::cout << "random seed = " << pars.meta.randomSeed + currentSlice * pars.numPlanes << std::endl;
					srand(pars.meta.randomSeed + currentSlice * pars.numPlanes);
					std::default_random_engine de(pars.meta.randomSeed + currentSlice * pars.numPlanes);
					normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0, 1);

					for (auto atom_num = 0; atom_num < x.size(); ++atom_num)
					{
						if (zPlane[atom_num] == currentSlice)
						{
							if (pars.meta.includeOccupancy)
							{
								if (static_cast<PRISMATIC_FLOAT_PRECISION>(rand()) / static_cast<PRISMATIC_FLOAT_PRECISION>(RAND_MAX) > occ[atom_num])
								{
									continue;
								}
							}
							//								if ( !pars.meta.includeOccupancy || static_cast<PRISMATIC_FLOAT_PRECISION>(rand())/static_cast<PRISMATIC_FLOAT_PRECISION> (RAND_MAX) <= occ[atom_num]) {
							const size_t cur_Z = Z_lookup[ID[atom_num]];
							PRISMATIC_FLOAT_PRECISION X, Y;
							if (pars.meta.includeThermalEffects)
							{ // apply random perturbations
								X = round((x[atom_num] + randn(de) * sigma[atom_num]) / pars.pixelSize[1]);
								Y = round((y[atom_num] + randn(de) * sigma[atom_num]) / pars.pixelSize[0]);
							}
							else
							{
								X = round((x[atom_num]) / pars.pixelSize[1]); // this line uses no thermal factor
								Y = round((y[atom_num]) / pars.pixelSize[0]); // this line uses no thermal factor
							}
							xp = xvec + (long)X;
							for (auto &i : xp)
								i = (i % dim1 + dim1) % dim1; // make sure to get a positive value

							yp = yvec + (long)Y;
							for (auto &i : yp)
								i = (i % dim0 + dim0) % dim0; // make sure to get a positive value
							for (auto ii = 0; ii < xp.size(); ++ii)
							{
								for (auto jj = 0; jj < yp.size(); ++jj)
								{
									// fill in value with lookup table
									projectedPotential.at(yp[jj], xp[ii]) += potentialLookup.at(cur_Z, jj, ii);
								}
							}
							//								}
						}
					}
					// copy the result to the full array
					copy(projectedPotential.begin(), projectedPotential.end(), &pars.pot.at(currentSlice, 0, 0));
					#ifdef PRISMATIC_BUILDING_GUI
					pars.progressbar->signalPotentialUpdate(currentSlice, pars.numPlanes);
					#endif //PRISMATIC_BUILDING_GUI
					++currentSlice;
				}
			}
		}));
	}
	cout << "Waiting for threads...\n";
	for (auto &t : workers)
		t.join();
#ifdef PRISMATIC_BUILDING_GUI
	pars.progressbar->setProgress(100);
#endif //PRISMATIC_BUILDING_GUI
};

void interpolatePotential(Array3D<PRISMATIC_FLOAT_PRECISION> &potShift,
							const Array3D<PRISMATIC_FLOAT_PRECISION> &potCrop,
							const PRISMATIC_FLOAT_PRECISION &wx,
							const PRISMATIC_FLOAT_PRECISION &wy,
							const PRISMATIC_FLOAT_PRECISION &wz,
							const size_t &xind,
							const size_t &yind,
							const size_t &zind)
{
	for(auto k = 0; k < potCrop.get_dimk(); k++)
	{
		for(auto j = 0; j < potCrop.get_dimj(); j++)
		{
			for(auto i = 0; i < potCrop.get_dimj(); i++)
			{
				potShift.at(k+zind,j+yind,i+xind) += potCrop.at(k,j,i)*wx*wy*wz;
			}
		}
	}
};

void cropLookup(Array3D<PRISMATIC_FLOAT_PRECISION> &potCrop,
				const Array4D<PRISMATIC_FLOAT_PRECISION> &potLookup,
				const size_t &cur_Z)
{
	//crops faces off of potLookup
	for(auto k = 0; k < potCrop.get_dimk(); k++)
	{
		for(auto j = 0; j < potCrop.get_dimj(); j++)
		{
			for(auto i = 0; i < potCrop.get_dimi(); i++)
			{
				potCrop.at(k,j,i) = potLookup.at(cur_Z, k+1, j+1, i+1);
			}
		}
	}

};			

void generateProjectedPotentials3D(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
								   const Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> &potLookup,
								   const vector<size_t> &unique_species,
								   const Array1D<long> &xvec,
								   const Array1D<long> &yvec,
								   const Array1D<PRISMATIC_FLOAT_PRECISION> &zvec)
{		
	long numPlanes = round(pars.tiledCellDim[0]/pars.meta.sliceThickness);
	//check if intermediate output was specified, if so, create index of output slices
	pars.numPlanes = numPlanes;
	if (pars.meta.numSlices == 0) pars.numSlices = pars.numPlanes;

	pars.pot = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{ (size_t) numPlanes, pars.imageSize[0], pars.imageSize[1]}});

	// create arrays for the coordinates
	Array1D<PRISMATIC_FLOAT_PRECISION> x = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> y = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> z = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> ID = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> sigma = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> occ = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.atoms.size()}});


	// populate arrays from the atoms structure
	for (auto i = 0; i < pars.atoms.size(); ++i)
	{
		x[i] = pars.atoms[i].x * pars.tiledCellDim[2];
		y[i] = pars.atoms[i].y * pars.tiledCellDim[1];
		z[i] = pars.atoms[i].z * pars.tiledCellDim[0];
		ID[i] = pars.atoms[i].species;
		sigma[i] = pars.atoms[i].sigma;
		occ[i] = pars.atoms[i].occ;
	}

	const long dim1 = (long) pars.pot.get_dimi();
	const long dim0 = (long) pars.pot.get_dimj();

	// correct z orientation
	auto max_z = pars.tiledCellDim[0];

	std::transform(z.begin(), z.end(), z.begin(), [&max_z](PRISMATIC_FLOAT_PRECISION &t_z) {
		return (-t_z + max_z); // If the +0.5 was to make the first slice z=1 not 0, can drop the +0.5 and -1
	});

	Array1D<PRISMATIC_FLOAT_PRECISION> zr = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{zvec.get_dimi()}});
	for (auto j = 0; j < zr.size(); ++j) zr[j] = (PRISMATIC_FLOAT_PRECISION)zvec[j] * pars.dzPot;

	//create fourier coordinate arrays for use in fourier shifting
	const PRISMATIC_FLOAT_PRECISION pi = std::acos(-1);
	Array1D<PRISMATIC_FLOAT_PRECISION> qy = makeFourierCoords(yvec.size(), (PRISMATIC_FLOAT_PRECISION) 1.0);
	Array1D<PRISMATIC_FLOAT_PRECISION> qx = makeFourierCoords(xvec.size(), (PRISMATIC_FLOAT_PRECISION) 1.0);
	std::pair<Array2D<PRISMATIC_FLOAT_PRECISION>, Array2D<PRISMATIC_FLOAT_PRECISION>> qmesh = meshgrid(qy,qx);
	Array2D<PRISMATIC_FLOAT_PRECISION> qya = qmesh.first;
	Array2D<PRISMATIC_FLOAT_PRECISION> qxa = qmesh.second;


	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> qyShift = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{qya.get_dimj(), qya.get_dimi()}});
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> qxShift = zeros_ND<2,std::complex<PRISMATIC_FLOAT_PRECISION>>({{qya.get_dimj(), qya.get_dimi()}});
	Array2D<PRISMATIC_FLOAT_PRECISION> q1(qya);
	for(auto j = 0; j < qya.get_dimj(); j++)
	{
		for(auto i = 0; i < qya.get_dimi(); i++)
		{
			qyShift.at(j,i) = -2.0*i*pi*qya.at(j,i);
			qxShift.at(j,i) = -2.0*i*pi*qxa.at(j,i);
		}
	}

	//band limit arrays for soft aperture with shift in realpace and fourier space; hard-coded for now
	PRISMATIC_FLOAT_PRECISION rband_max = 0.95;
	PRISMATIC_FLOAT_PRECISION rband_min = 0.75;
	Array2D<PRISMATIC_FLOAT_PRECISION> qband(q1);
	qband *= -2.0;
	qband += rband_max;
	qband /= (rband_max-rband_min);
	for(auto i = 0; i < qband.size(); i++)
	{
		qband[i] = (rband_max-2*q1[i])/(rband_max-rband_min);
		qband[i] = std::max(qband[i], (PRISMATIC_FLOAT_PRECISION) 0.0);
		qband[i] = std::min(qband[i], (PRISMATIC_FLOAT_PRECISION) 1.0);
		qband[i] = pow(std::sin(qband[i]*pi/2.0), 2.0);
	}

	std::pair<Array2D<long>, Array2D<long>> rmesh = meshgrid(yvec,xvec);
	Array2D<PRISMATIC_FLOAT_PRECISION> rband(qband); //construct with qband to avoid type mismatch, knowing sizes are the same
	PRISMATIC_FLOAT_PRECISION xl = (PRISMATIC_FLOAT_PRECISION) xvec[xvec.size()-1];
	PRISMATIC_FLOAT_PRECISION yl = (PRISMATIC_FLOAT_PRECISION) yvec[yvec.size()-1];
	for(auto i =0; i < rband.size(); i++)
	{
		rband[i] = pow((rmesh.first[i] / (yl+0.5)), 2.0) + pow((rmesh.second[i] / (xl+0.5)), 2.0);
		rband[i] = (rband[i] <= 1) ? 1.0 : 0.0;  
	}


	// create a key-value map to match the atomic Z numbers with their place in the potential lookup table
	map<size_t, size_t> Z_lookup;
	for (auto i = 0; i < unique_species.size(); ++i)
		Z_lookup[unique_species[i]] = i;
		
	std::vector<std::thread> workers;
	size_t numWorkers = pars.meta.numThreads; //std::min(pars.meta.numThreads, (size_t) 4); //heuristic for now, TODO: improve parallelization scheme to segment atoms over regions to avoid write locks
	workers.reserve(numWorkers);
	WorkDispatcher dispatcher(0, pars.atoms.size());
	const size_t print_frequency = std::max((size_t)1, pars.atoms.size() / 10);

	PRISMATIC_FFTW_INIT_THREADS();
	std::cout << "Base random seed = " << pars.meta.randomSeed << std::endl;
	for (long t = 0; t < numWorkers; t++)
	{
		std::cout << "Launching thread #" << t << " to compute projected potential slices\n";
		workers.push_back(thread([&pars, &x, &y, &z, &ID, &sigma, &occ, &print_frequency,
								 &Z_lookup, &xvec, &yvec, &zvec, &zr, &dim0, &dim1,
								 &numPlanes, &potLookup, &rband, &qband, &qxShift, &qyShift, &dispatcher]()
		{
			size_t currentAtom, stop;
			currentAtom = stop = 0;
			while (dispatcher.getWork(currentAtom, stop))
			{
				while(currentAtom != stop)
				{
					if(!(currentAtom % print_frequency))
					{
						std::cout << "Computing atom " << currentAtom << "/" << pars.atoms.size() << std::endl;
					}

					// create a random number generator to simulate thermal effects
					srand(pars.meta.randomSeed+currentAtom);
					std::default_random_engine de(pars.meta.randomSeed+currentAtom);
					normal_distribution<PRISMATIC_FLOAT_PRECISION> randn(0, 1);
					
					const size_t cur_Z = Z_lookup[ID[currentAtom]];
					PRISMATIC_FLOAT_PRECISION X, Y, Z;
					PRISMATIC_FLOAT_PRECISION perturbX, perturbY, perturbZ;
					if (pars.meta.includeThermalEffects)
					{ // apply random perturbations
						perturbX = randn(de) * sigma[currentAtom];
						perturbY = randn(de) * sigma[currentAtom];
						perturbZ = randn(de) * sigma[currentAtom];
						X = round((x[currentAtom] + perturbX) / pars.pixelSize[1]);
						Y = round((y[currentAtom] + perturbY) / pars.pixelSize[0]);
						Z = (z[currentAtom] + perturbZ); //z gets rounded and normalized later
					}
					else
					{
						perturbX = perturbY = perturbZ = 0;
						X = round((x[currentAtom]) / pars.pixelSize[1]); // this line uses no thermal factor
						Y = round((y[currentAtom]) / pars.pixelSize[0]); // this line uses no thermal factor
						Z = (z[currentAtom]); // this line uses no thermal factor, z gets rounded and normalized later
					}

					PRISMATIC_FLOAT_PRECISION dxPx = (x[currentAtom] + perturbX)/ pars.pixelSize[1] - X;
					PRISMATIC_FLOAT_PRECISION dyPy = (y[currentAtom] + perturbY)/ pars.pixelSize[0] - Y;

					Array1D<long> xp = xvec + (long) X;
					Array1D<long> yp = yvec + (long) Y;

					for(auto &i : xp) i = (i % dim1 + dim1) % dim1;
					for(auto &i : yp) i = (i % dim0 + dim0) % dim0;
					Array1D<long> zp = zeros_ND<1, long>({{zvec.get_dimi()}});
					std::vector<long> zVals(zp.size(), 0);
					for(auto i = 0; i < zp.size(); i++)
					{
						PRISMATIC_FLOAT_PRECISION tmp = round((Z+zr[i])/pars.meta.sliceThickness + 0.5)-1;
						tmp = std::max(tmp, (PRISMATIC_FLOAT_PRECISION) 0.0);
						zp[i] = std::min((long) tmp, numPlanes-1);
						zVals[i] = zp[i];
					}

					std::sort(zVals.begin(), zVals.end());
					auto last = std::unique(zVals.begin(), zVals.end());
					zVals.erase(last, zVals.end());

					//iterate through unique z slice values
					for(auto cz_ind = 0; cz_ind < zVals.size(); cz_ind++)
					{
						
						//create tmp array to add potential lookup table to
						Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> tmp_pot = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{yp.size(), xp.size()}});

						for(auto kk = 0; kk < zp.size(); kk++)
						{
							if(zp[kk] == zVals[cz_ind])
							{
								for(auto jj = 0; jj < yp.size(); jj++)
								{
									for(auto ii = 0; ii < xp.size(); ii++)
									{
										tmp_pot.at(jj,ii) += potLookup.at(cur_Z, kk,jj,ii);
									}
								}
							}
						}

						//apply fourier shift and qband limit
						for(auto jj = 0; jj < yp.size(); jj++)
						{
							for(auto ii = 0; ii < xp.size(); ii++)
							{
								tmp_pot.at(jj,ii) *= qband.at(jj,ii) * exp(qxShift.at(jj,ii)*dxPx + qyShift.at(jj,ii)*dyPy);
							}
						}

						//inverse FFT and normalize by size of array
						unique_lock<mutex> gatekeeper(fftw_plan_lock);
						PRISMATIC_FFTW_PLAN plan_inverse = PRISMATIC_FFTW_PLAN_DFT_2D(tmp_pot.get_dimj(), tmp_pot.get_dimi(),
																				reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&tmp_pot[0]),
																				reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&tmp_pot[0]),
																				FFTW_BACKWARD,
																				FFTW_ESTIMATE);
						gatekeeper.unlock();
						PRISMATIC_FFTW_EXECUTE(plan_inverse);
						gatekeeper.lock();
						PRISMATIC_FFTW_DESTROY_PLAN(plan_inverse);
						gatekeeper.unlock();
						for(auto &t : tmp_pot) t /= tmp_pot.get_dimi()*tmp_pot.get_dimj();

						//apply realspace band limit
						for(auto i = 0; i < tmp_pot.size(); i++) tmp_pot[i] *= rband[i];

						//then write
						//put into a mutex lock to prevent race condition on potential writing when atoms overlap within potential bound
						std::unique_lock<std::mutex> write_gatekeeper(potentialWriteLock);
						for(auto jj = 0; jj < yp.size(); jj++)
						{
							for(auto ii = 0; ii < xp.size(); ii++)
							{
								pars.pot.at(zVals[cz_ind],yp[jj],xp[ii]) += tmp_pot.at(jj,ii).real();
							}
						}
						write_gatekeeper.unlock();
					}
					++currentAtom;
				}
			}
		}));
	}
	std::cout << "Waiting for threads...\n";
	for (auto &t : workers)
		t.join();

	PRISMATIC_FFTW_CLEANUP_THREADS();
};

void PRISM01_calcPotential(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	//builds projected, sliced potential
	
	// setup some coordinates
	cout << "Entering PRISM01_calcPotential" << endl;
	PRISMATIC_FLOAT_PRECISION yleng = std::ceil(pars.meta.potBound / pars.pixelSize[0]);
	PRISMATIC_FLOAT_PRECISION xleng = std::ceil(pars.meta.potBound / pars.pixelSize[1]);
	ArrayND<1, vector<long>> xvec(vector<long>(2 * (size_t)xleng + 1, 0), {{2 * (size_t)xleng + 1}});
	ArrayND<1, vector<long>> yvec(vector<long>(2 * (size_t)yleng + 1, 0), {{2 * (size_t)yleng + 1}});
	{
		PRISMATIC_FLOAT_PRECISION tmpx = -xleng;
		PRISMATIC_FLOAT_PRECISION tmpy = -yleng;
		for (auto &i : xvec)
			i = tmpx++;
		for (auto &j : yvec)
			j = tmpy++;
	}
	Array1D<PRISMATIC_FLOAT_PRECISION> xr(vector<PRISMATIC_FLOAT_PRECISION>(2 * (size_t)xleng + 1, 0), {{2 * (size_t)xleng + 1}});
	Array1D<PRISMATIC_FLOAT_PRECISION> yr(vector<PRISMATIC_FLOAT_PRECISION>(2 * (size_t)yleng + 1, 0), {{2 * (size_t)yleng + 1}});
	for (auto i = 0; i < xr.size(); ++i)
		xr[i] = (PRISMATIC_FLOAT_PRECISION)xvec[i] * pars.pixelSize[1];
	for (auto j = 0; j < yr.size(); ++j)
		yr[j] = (PRISMATIC_FLOAT_PRECISION)yvec[j] * pars.pixelSize[0];

	vector<size_t> unique_species = get_unique_atomic_species(pars);

	if(pars.meta.potential3D)
	{	//set up Z coords

		pars.dzPot = pars.meta.sliceThickness/pars.meta.zSampling;
        PRISMATIC_FLOAT_PRECISION zleng = std::ceil(pars.meta.potBound/pars.dzPot);
		Array1D<PRISMATIC_FLOAT_PRECISION> zvec = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{(size_t)zleng*2}});
		for (auto j = -zleng; j < zleng; j++)
		{
			zvec[j+zleng] = (PRISMATIC_FLOAT_PRECISION) j + 0.5;
		}

		Array1D<PRISMATIC_FLOAT_PRECISION> zr(zvec);
        for (auto j = 0; j < zr.size(); ++j) zr[j] = zvec[j] * pars.dzPot;

		// initialize the lookup table and precompute unique potentials
		Array4D<std::complex<PRISMATIC_FLOAT_PRECISION>> potentialLookup = zeros_ND<4, std::complex<PRISMATIC_FLOAT_PRECISION>>({{unique_species.size(), 2 * (size_t)zleng, 2 * (size_t)yleng + 1, 2 * (size_t)xleng + 1}});
		fetch_potentials3D(potentialLookup, unique_species, xr, yr, zr);
		//generate potential
		generateProjectedPotentials3D(pars, potentialLookup, unique_species, xvec, yvec, zvec);

	}else{
		// initialize the lookup table
		Array3D<PRISMATIC_FLOAT_PRECISION> potentialLookup = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{unique_species.size(), 2 * (size_t)yleng + 1, 2 * (size_t)xleng + 1}});

		// precompute the unique potentials
		fetch_potentials(potentialLookup, unique_species, xr, yr);

		// populate the slices with the projected potentials
		generateProjectedPotentials(pars, potentialLookup, unique_species, xvec, yvec);
	}

	if (pars.meta.savePotentialSlices) 
	{
		std::cout << "Writing potential slices to output file." << std::endl;
		savePotentialSlices(pars);
	}
}

void PRISM01_importPotential(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	std::cout << "Setting up PRISM01 auxilary variables according to " << pars.meta.importFile << " metadata." << std::endl;
	//scope out imported tmp_pot as soon as possible
	{
		Array3D<PRISMATIC_FLOAT_PRECISION> tmp_pot;
		if(pars.meta.importPath.size() > 0)
		{
			readRealDataSet_inOrder(tmp_pot, pars.meta.importFile, pars.meta.importPath);
		}
		else //read default path
		{
			std::string groupPath = "4DSTEM_simulation/data/realslices/ppotential_fp" + getDigitString(pars.fpFlag) + "/data";
			readRealDataSet_inOrder(tmp_pot, pars.meta.importFile, groupPath);
		}

		//initailize array and get data in right order
		pars.pot = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{tmp_pot.get_dimi(), tmp_pot.get_dimj(), tmp_pot.get_dimk()}});
		for(auto i = 0; i < tmp_pot.get_dimi(); i++)
		{
			for(auto j = 0; j < tmp_pot.get_dimj(); j++)
			{
				for(auto k = 0; k < tmp_pot.get_dimk(); k++)
				{
					pars.pot.at(i,j,k) = tmp_pot.at(k,j,i);
				}
			}
		}
	}

	pars.numPlanes = pars.pot.get_dimk();
	if (pars.meta.numSlices == 0)
	{
		pars.numSlices = pars.numPlanes;
	}

	//resample coordinates if PRISM algorithm and size of PS array in not a multiple of 4*fx or 4*fy
	if(pars.meta.algorithm == Algorithm::PRISM)
	{
		if ( (pars.pot.get_dimi() % 4*pars.meta.interpolationFactorX) || (pars.pot.get_dimj() % pars.meta.interpolationFactorY))
		{
			std::cout << "Resampling imported potential to align grid size with requested interpolation factors fx = " 
					  << pars.meta.interpolationFactorX << " and fy = " << pars.meta.interpolationFactorY << std::endl;
			fourierResampling(pars);
		}
	}

	//TODO: metadata from non-prismatic sources?
    std::string groupPath = "4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters";
	PRISMATIC_FLOAT_PRECISION meta_cellDims[3];
	readAttribute(pars.meta.importFile, groupPath, "c", meta_cellDims);

	PRISMATIC_FLOAT_PRECISION meta_tile[3];
	readAttribute(pars.meta.importFile, groupPath, "t", meta_tile);

	pars.tiledCellDim[0] = meta_cellDims[2]*meta_tile[2];
	pars.tiledCellDim[1] = meta_cellDims[1]*meta_tile[1];
	pars.tiledCellDim[2] = meta_cellDims[0]*meta_tile[0];

	std::vector<PRISMATIC_FLOAT_PRECISION> pixelSize{(PRISMATIC_FLOAT_PRECISION) pars.tiledCellDim[1], (PRISMATIC_FLOAT_PRECISION) pars.tiledCellDim[2]};
	pars.imageSize[0] = pars.pot.get_dimj();
	pars.imageSize[1] = pars.pot.get_dimi();
	pixelSize[0] /= pars.imageSize[0];
	pixelSize[1] /= pars.imageSize[1];
	pars.pixelSize = pixelSize;

	if (pars.meta.savePotentialSlices) 
	{
		std::cout << "Writing potential slices to output file." << std::endl;
		savePotentialSlices(pars);
	}

};

void fourierResampling(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	int Ni = 0;
	int Nj = 0;

	//get highest multiple of 4*fx and 4*fy to ensure resampling to a smaller grid only
	while(Ni < pars.pot.get_dimi()) Ni += pars.meta.interpolationFactorX*4;
	while(Nj < pars.pot.get_dimj()) Nj += pars.meta.interpolationFactorY*4;
	Ni -= pars.meta.interpolationFactorX*4;
 	Nj -= pars.meta.interpolationFactorY*4;
 	
	Array3D<PRISMATIC_FLOAT_PRECISION> newPot = zeros_ND<3,PRISMATIC_FLOAT_PRECISION>({{pars.pot.get_dimk(), (size_t) Nj, (size_t) Ni}});

	//create storage variables to hold data from FFTs
	Array2D<complex<PRISMATIC_FLOAT_PRECISION>> fstore = zeros_ND<2,complex<PRISMATIC_FLOAT_PRECISION>>({{pars.pot.get_dimj(), pars.pot.get_dimi()}});
	Array2D<complex<PRISMATIC_FLOAT_PRECISION>> bstore = zeros_ND<2,complex<PRISMATIC_FLOAT_PRECISION>>({{(size_t) Nj, (size_t) Ni}});
	Array2D<complex<PRISMATIC_FLOAT_PRECISION>> fpot = zeros_ND<2,complex<PRISMATIC_FLOAT_PRECISION>>({{pars.pot.get_dimj(),pars.pot.get_dimi()}});
	Array2D<complex<PRISMATIC_FLOAT_PRECISION>> bpot = zeros_ND<2,complex<PRISMATIC_FLOAT_PRECISION>>({{(size_t)Nj,(size_t) Ni}});
	
	//create FFT plans 
	PRISMATIC_FFTW_INIT_THREADS();
	PRISMATIC_FFTW_PLAN_WITH_NTHREADS(pars.meta.numThreads);
	
	unique_lock<mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan_forward = PRISMATIC_FFTW_PLAN_DFT_2D(fstore.get_dimj(), fstore.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&fpot[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&fstore[0]),
															FFTW_FORWARD,
															FFTW_ESTIMATE);

	PRISMATIC_FFTW_PLAN plan_inverse = PRISMATIC_FFTW_PLAN_DFT_2D(bstore.get_dimj(), bstore.get_dimi(),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&bstore[0]),
															reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&bpot[0]),
															FFTW_BACKWARD,
															FFTW_ESTIMATE);
	gatekeeper.unlock();

	//calculate indices for downsampling in fourier space
	int nyqi = std::floor(Ni/2) + 1;
	int nyqj = std::floor(Nj/2) + 1;

	for(auto k = 0; k < newPot.get_dimk(); k++)
	{
		//copy current slice to forward transform
		for(auto i = 0; i < fpot.size(); i++) fpot[i] = pars.pot[k*pars.pot.get_dimj()*pars.pot.get_dimi()+i];
		
		//forward transform 
		PRISMATIC_FFTW_EXECUTE(plan_forward);

		//copy relevant quadrants to backward store
		//manual looping through quadrants
		for(auto j = 0; j < nyqj; j++)
		{
			for(auto i = 0; i < nyqi; i++)
			{
				bstore.at(j, i) = fstore.at(j, i);
			}
		}

		for(auto j = nyqj-Nj; j < 0; j++)
		{
			for(auto i = 0; i < nyqi; i++)
			{
				bstore.at(Nj + j, i) = fstore.at(fstore.get_dimj() + j, i);
			}
		}

		for(auto j = 0; j < nyqj; j++)
		{
			for(auto i = nyqi-Ni; i < 0; i++)
			{
				bstore.at(j, Ni + i) = fstore.at(j, fstore.get_dimi() + i);
			}
		}

		for(auto j = nyqj-Nj; j < 0; j++)
		{
			for(auto i = nyqi-Ni; i < 0; i++)
			{
				bstore.at(Nj + j, Ni + i) = fstore.at(fstore.get_dimj() + j, fstore.get_dimi() + i);
			}
		}

		//inverse transform
		PRISMATIC_FFTW_EXECUTE(plan_inverse);

		//store slice in potential
		for(auto i = 0; i < bpot.size(); i++) newPot[k*newPot.get_dimj()*newPot.get_dimi()+i] = bpot[i].real();
	}

	//store final resort after normalizing FFT, rescaling from transform, and removing negative values
	PRISMATIC_FLOAT_PRECISION orig_x = pars.pot.get_dimi();
	PRISMATIC_FLOAT_PRECISION orig_y = pars.pot.get_dimj();
	PRISMATIC_FLOAT_PRECISION new_x = Ni;
	PRISMATIC_FLOAT_PRECISION new_y = Nj;
	newPot /= Ni*Nj;
	newPot *= (new_x/orig_x)*(new_y/orig_y);

	pars.pot = newPot;
};

} // namespace Prismatic
