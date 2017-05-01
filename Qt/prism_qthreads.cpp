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
//    this->parent->potential = params.pot;
    // indicate that the potential is ready
    this->parent->potentialReady = true;
    if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");

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
//        QMutexLocker gatekeeper(&this->parent->dataLock);

        // perform copy
//        this->parent->pars = params;
//        this->parent->potential = params.pot;
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


    std::cout << "calculating S-Matrix" << std::endl;
//    // calculate S-Matrix
    PRISM::PRISM02(params);
////    QMutexLocker gatekeeper(&this->parent->sMatrixLock);
    QMutexLocker gatekeeper(&this->parent->dataLock);

//    // perform copy
    this->parent->pars = params;
    this->parent->ScompactReady = true;
}


ProbeThread::ProbeThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}

void ProbeThread::run(){};

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

    std::cout << "Calculation complete" << std::endl;
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


    std::cout << "Calculation complete" << std::endl;
    std::cout<<"after copy this->parent->output.at(0,0,0) = " << this->parent->pars.output.at(0,0,0) << std::endl;

}

PotentialThread::~PotentialThread(){}
SMatrixThread::~SMatrixThread(){}
ProbeThread::~ProbeThread(){}
FullPRISMCalcThread::~FullPRISMCalcThread(){}
FullMultisliceCalcThread::~FullMultisliceCalcThread(){}
