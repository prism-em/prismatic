// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "prism_qthreads.h"
#include "prism_progressbar.h"
#include "PRISM01.h"
#include "PRISM02.h"
#include "PRISM03.h"
#include "Multislice.h"
#include "utility.h"
#include <iostream>

PotentialThread::PotentialThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}

void PotentialThread::run(){
    // create parameters
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
    QMutexLocker calculationLocker(&this->parent->calculationLock);
    // calculate potential
    PRISM::PRISM01(params);
    // acquire the mutex so we can safely copy to the GUI copy of the potential
    QMutexLocker gatekeeper(&this->parent->dataLock);
    // perform copy
    this->parent->pars = params;
    // indicate that the potential is ready
    this->parent->potentialReady = true;
    if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
    std::cout << "Projected potential calculation complete" << std::endl;
}

SMatrixThread::SMatrixThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
        parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}

void SMatrixThread::run(){
    QMutexLocker calculationLocker(&this->parent->calculationLock);
    PRISM::configure(meta);
    // create parameters
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
    // calculate potential if it hasn't been already
    if (!this->parent->potentialReady){
        // calculate potential
        PRISM::PRISM01(params);
        // indicate that the potential is ready
        this->parent->potentialReady = true;
        if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
        emit potentialCalculated();
    } else {
        std::cout << "COPYING" << std::endl;
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "Potential already calculated. Using existing result." << std::endl;
    }


//    // calculate S-Matrix
    PRISM::PRISM02(params);
////    QMutexLocker gatekeeper(&this->parent->sMatrixLock);
    QMutexLocker gatekeeper(&this->parent->dataLock);

//    // perform copy
    this->parent->pars = params;
    this->parent->ScompactReady = true;
    std::cout << "S-matrix calculation complete" << std::endl;
}


ProbeThread::ProbeThread(PRISMMainWindow *_parent, PRISM_FLOAT_PRECISION _X, PRISM_FLOAT_PRECISION _Y, prism_progressbar *_progressbar) :
parent(_parent), X(_X), Y(_Y), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}

void ProbeThread::run(){
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
//    PRISM::Parameters<PRISM_FLOAT_PRECISION> params_multi(params);

    QMutexLocker calculationLocker(&this->parent->calculationLock);

    PRISM::configure(meta);
    if (!this->parent->potentialReady){
    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
//        params_multi = params;
        this->parent->pars = params;
//        this->parent->pars_multi = params_multi;
        this->parent->potentialReady = true;
        if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
    }
    emit potentialCalculated();
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
//        params_multi = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "Potential already calculated. Using existing result." << std::endl;
    }

    if (!this->parent->ScompactReady){
    PRISM::PRISM02(params);
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);

        // perform copy
        this->parent->pars = params;
//        this->parent->pars_multi = params_multi;
        // indicate that the potential is ready
        this->parent->ScompactReady = true;
    }
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
//        params_multi = this->parent->pars_multi;
        params.progressbar = progressbar;
        std::cout << "S-Matrix already calculated. Using existing result." << std::endl;
    }

    std::pair<PRISM::Array2D< std::complex<PRISM_FLOAT_PRECISION> >, PRISM::Array2D< std::complex<PRISM_FLOAT_PRECISION> > > prism_probes, multislice_probes;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params_multi(params);

    // setup and calculate PRISM probe
    PRISM::setupCoordinates_2(params);

    // setup angles of detector and image sizes
    PRISM::setupDetector(params);

    // setup coordinates and indices for the beams
    PRISM::setupBeams_2(params);

    // setup Fourier coordinates for the S-matrix
    PRISM::setupFourierCoordinates(params);

    // initialize the output to the correct size for the output mode
    PRISM::createStack_integrate(params);

//  perform some necessary setup transformations of the data
    PRISM::transformIndices(params);

    // initialize/compute the probes
    PRISM::initializeProbes(params);

    prism_probes = PRISM::getSinglePRISMProbe_CPU(params, X, Y);


    // setup and calculate Multislice probe
    // setup coordinates and build propagators
    PRISM::setupCoordinates_multislice(params_multi);

    // setup detector coordinates and angles
    PRISM::setupDetector_multislice(params_multi);

    // create initial probes
    PRISM::setupProbes_multislice(params_multi);

    // create transmission array
    PRISM::createTransmission(params_multi);

    // initialize output stack
    PRISM::createStack(params_multi);

    multislice_probes = PRISM::getSingleMultisliceProbe_CPU(params_multi, X, Y);

    PRISM::Array2D<PRISM_FLOAT_PRECISION> debug = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    for (auto j = 0; j < multislice_probes.first.get_dimj(); ++j){
        for (auto i = 0; i < multislice_probes.first.get_dimi(); ++i){
        debug.at(j,i) = std::abs(multislice_probes.first.at(j,i));
        }
    }
    debug.toMRC_f("/mnt/spareA/clion/PRISM/build/db.mrc");
    for (auto j = 0; j < multislice_probes.second.get_dimj(); ++j){
        for (auto i = 0; i < multislice_probes.second.get_dimi(); ++i){
        debug.at(j,i) = std::abs(multislice_probes.second.at(j,i));
        }
    }
    debug.toMRC_f("/mnt/spareA/clion/PRISM/build/dbk.mrc");

    QMutexLocker gatekeeper(&this->parent->dataLock);
    // perform copy
    this->parent->pars = params;
    this->parent->probeSetupReady = true;

    prism_probes = upsamplePRISMProbe(prism_probes.first, multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi());

//    for (auto j = 0; j < multislice_probes.first.get_dimj(); ++j){
//        for (auto i = 0; i < multislice_probes.first.get_dimi(); ++i){
//        debug.at(j,i) = std::abs(multislice_probes.first.at(j,i));
//        }
//    }
//    debug.toMRC_f("/mnt/spareA/clion/PRISM/build/db.mrc");
//    for (auto j = 0; j < multislice_probes.second.get_dimj(); ++j){
//        for (auto i = 0; i < multislice_probes.second.get_dimi(); ++i){
//        debug.at(j,i) = std::abs(multislice_probes.second.at(j,i));
//        }
//    }
//    debug.toMRC_f("/mnt/spareA/clion/PRISM/build/dbk.mrc");


    std::cout << "prism_probes.first.get_dimj() = " << prism_probes.first.get_dimj() <<std::endl;
    std::cout << "multislice_probes.first.get_dimj() = " << multislice_probes.first.get_dimj() <<std::endl;

    PRISM::Array2D<PRISM_FLOAT_PRECISION> pr    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{prism_probes.first.get_dimj(), prism_probes.first.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> pk    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{prism_probes.second.get_dimj(), prism_probes.second.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> mr    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> mk    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.second.get_dimj(), multislice_probes.second.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> diffr = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> diffk = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.second.get_dimj(), multislice_probes.second.get_dimi()}});

//if (parent->ui->checkBox_log->checked()){
//    for (auto i = 0; i < prism_probes.first.size(); ++i){
//        pr[i] =  std::log(std::abs(prism_probes.first[i]));
//    }
//    for (auto i = 0; i < prism_probes.second.size(); ++i){
//        pk[i] =  std::log(std::abs(prism_probes.second[i]));
//    }
//    for (auto i = 0; i < multislice_probes.first.size(); ++i){
//        mr[i] =  std::log(std::abs(multislice_probes.first[i]));
//    }
//    for (auto i = 0; i < multislice_probes.second.size(); ++i){
//        mk[i] =  std::log(std::abs(multislice_probes.second[i]));
//    }
//} else{
    for (auto i = 0; i < prism_probes.first.size(); ++i){
        pr[i] =  std::abs(prism_probes.first[i]);
    }
    for (auto i = 0; i < prism_probes.second.size(); ++i){
        pk[i] = std::abs(prism_probes.second[i]);
    }
    for (auto i = 0; i < multislice_probes.first.size(); ++i){
        mr[i] =  std::abs(multislice_probes.first[i]);
    }
    for (auto i = 0; i < multislice_probes.second.size(); ++i){
        mk[i] =  std::abs(multislice_probes.second[i]);
    }
//}

    for (auto i = 0; i < prism_probes.second.size(); ++i){
        diffr[i] =  (std::abs(pr[i] - mr[i]));
        diffk[i] =  (std::abs(pk[i] - mk[i]));
    }

//    for (auto i = 0; i < prism_probes.first.size(); ++i){
//        pr[i] =  std::abs(prism_probes.first[i]);
//    }
//    for (auto i = 0; i < prism_probes.second.size(); ++i){
//        pk[i] =  std::abs(prism_probes.second[i]);
//    }
//    for (auto i = 0; i < multislice_probes.first.size(); ++i){
//        mr[i] =  std::abs(multislice_probes.first[i]);
//    }prism_probes
//    for (auto i = 0; i < multislice_probes.second.size(); ++i){
//        mk[i] =  std::abs(multislice_probes.second[i]);
//    }
//    for (auto i = 0; i < prism_probes.second.size(); ++i){
//        diffr[i] =  std::abs(pr[i] - mr[i]);
//        diffk[i] =  std::abs(pk[i] - mk[i]);
//    }
//    emit signalProbeR_PRISM(pr);
//    emit signalProbeK_PRISM(pk);
//    emit signalProbeR_Multislice(mr);
//    emit signalProbeK_Multislice(mk);
//    emit signalProbe_diffR(diffr);
//    emit signalProbe_diffK(diffk);
    emit signalProbeR_PRISM((pr));
    emit signalProbeK_PRISM(fftshift2(pk));
    emit signalProbeR_Multislice((mr));
    emit signalProbeK_Multislice(fftshift2(mk));
    emit signalProbe_diffR((diffr));
    emit signalProbe_diffK(fftshift2(diffk));
}

FullPRISMCalcThread::FullPRISMCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}


void FullPRISMCalcThread::run(){
    std::cout << "Full PRISM Calculation thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
    QMutexLocker calculationLocker(&this->parent->calculationLock);

    PRISM::configure(meta);
//  //  PRISM::Parameters<PRISM_FLOAT_PRECISION> params = PRISM::execute_plan(meta);
    if (!this->parent->potentialReady){
    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
    {
//        QMutexLocker gatekeeper(&this->parent->potentialLock);
        QMutexLocker gatekeeper(&this->parent->dataLock);
        this->parent->pars = params;
//        this->parent->potential = params.pot;
        this->parent->potentialReady = true;
        if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
    }
    emit potentialCalculated();
    } else {
        std::cout << "COPYING" << std::endl;
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "Potential already calculated. Using existing result." << std::endl;
    }

    if (!this->parent->ScompactReady){
    PRISM::PRISM02(params);
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);

        // perform copy
        this->parent->pars = params;
        // indicate that the potential is ready
        this->parent->ScompactReady = true;
    }
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "S-Matrix already calculated. Using existing result." << std::endl;
    }
//    emit ScompactCalculated();

    PRISM::PRISM03(params);
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        this->parent->pars = params;
        this->parent->detectorAngles = params.detectorAngles;
        for (auto& a:this->parent->detectorAngles) a*=1000; // convert to mrads
        this->parent->outputReady = true;
        params.output.toMRC_f(params.meta.filename_output.c_str());
    }
    emit outputCalculated();
    std::cout << "PRISM calculation complete" << std::endl;
}


FullMultisliceCalcThread::FullMultisliceCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
        parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}


void FullMultisliceCalcThread::run(){
    std::cout << "Full Multislice Calculation thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
	QMutexLocker calculationLocker(&this->parent->calculationLock);
    PRISM::configure(meta);
    if (!this->parent->potentialReady){
        PRISM::PRISM01(params);
        std::cout <<"Potential Calculated" << std::endl;
        {
    //        QMutexLocker gatekeeper(&this->parent->potentialLock);
            QMutexLocker gatekeeper(&this->parent->dataLock);
            this->parent->pars = params;
    //        this->parent->potential = params.pot;
            this->parent->potentialReady = true;
            if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
        }
        } else {
            std::cout << "COPYING" << std::endl;
            QMutexLocker gatekeeper(&this->parent->dataLock);
            params = this->parent->pars;
            params.progressbar = progressbar;
            std::cout << "Potential already calculated. Using existing result." << std::endl;
        }

    emit potentialCalculated();

    PRISM::Multislice(params);
    {
//        QMutexLocker gatekeeper(&this->parent->outputLock);
        QMutexLocker gatekeeper(&this->parent->dataLock);
        this->parent->pars = params;
	    this->parent->detectorAngles = params.detectorAngles;
	    for (auto& a:this->parent->detectorAngles) a*=1000; // convert to mrads
        this->parent->outputReady = true;
//        PRISM::Array3D<PRISM_FLOAT_PRECISION> reshaped_output = PRISM::zeros_ND<3, PRISM_FLOAT_PRECISION>(
//        {{params.output.get_diml(), params.output.get_dimk(), params.output.get_dimj()}});
//        auto ptr = reshaped_output.begin();
//        for (auto &i:params.output)*ptr++=i;
//        reshaped_output.toMRC_f(params.meta.filename_output.c_str());
        params.output.toMRC_f(params.meta.filename_output.c_str());
    }
    emit outputCalculated();
    std::cout << "Multislice calculation complete" << std::endl;
}

PotentialThread::~PotentialThread(){}
SMatrixThread::~SMatrixThread(){}
ProbeThread::~ProbeThread(){}
FullPRISMCalcThread::~FullPRISMCalcThread(){}
FullMultisliceCalcThread::~FullMultisliceCalcThread(){}
