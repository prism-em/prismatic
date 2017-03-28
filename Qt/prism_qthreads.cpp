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
    // calculate potential
    PRISM::PRISM01(params);
    // acquire the mutex so we can safely copy to the GUI copy of the potential
    QMutexLocker gatekeeper(&this->parent->potentialLock);
    // perform copy
    this->parent->potential = params.pot;
    // indicate that the potential is ready
    this->parent->potentialReady = true;
}

SMatrixThread::SMatrixThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
        parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}

void SMatrixThread::run(){
    // create parameters
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
    // calculate potential if it hasn't been already
    if (!this->parent->potentialReady){
        // calculate potential
        PRISM::PRISM01(params);
        // acquire the mutex so we can safely copy to the GUI copy of the potential
        QMutexLocker gatekeeper(&this->parent->potentialLock);
        // perform copy
        this->parent->potential = params.pot;
        // indicate that the potential is ready
        this->parent->potentialReady = true;
    }
    std::cout << "calculating S-Matrix" << std::endl;
    // calculate S-Matrix
    PRISM::PRISM02(params);
    // acquire the mutex so we can safely copy to the GUI copy of the potential
    QMutexLocker gatekeeper(&this->parent->sMatrixLock);
    // perform copy
    this->parent->Scompact = params.Scompact;
    // indicate that the potential is ready
    this->parent->ScompactReady = true;
    std::cout << "copying S-Matrix" << std::endl;
    std::cout << "S-Matrix.at(0,0,0) = " << this->parent->Scompact.at(0,0,0) << std::endl;
}


FullPRISMCalcThread::FullPRISMCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}


void FullPRISMCalcThread::run(){
    std::cout << "Full PRISM Calculation thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta);
    PRISM::configure(meta);
  //  PRISM::Parameters<PRISM_FLOAT_PRECISION> params = PRISM::execute_plan(meta);
    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
    {
        QMutexLocker gatekeeper(&this->parent->potentialLock);
        this->parent->potential = params.pot;
        this->parent->potentialReady = true;
    }


    PRISM::PRISM02(params);
    // acquire the mutex so we can safely copy to the GUI copy of the potential
    {
        QMutexLocker gatekeeper(&this->parent->sMatrixLock);
        // perform copy
        this->parent->Scompact = params.Scompact;
        // indicate that the potential is ready
        this->parent->ScompactReady = true;
    }
    std::cout << "copying S-Matrix" << std::endl;
    std::cout << "S-Matrix.at(0,0,0) = " << this->parent->Scompact.at(0,0,0) << std::endl;


    PRISM::PRISM03(params);
    {
        QMutexLocker gatekeeper(&this->parent->stackLock);
        this->parent->stack = params.stack;
        this->parent->stackReady = true;
    }

//    std::cout<<"after copy this->parent->pot.at(0,0,0) = " << this->parent->potential.at(0,0,0) << std::endl;
}


FullMultisliceCalcThread::FullMultisliceCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
        parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}


void FullMultisliceCalcThread::run(){
    std::cout << "Full Multislice Calculation thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta);
    PRISM::configure(meta);
    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
    {
        QMutexLocker gatekeeper(&this->parent->potentialLock);
        this->parent->potential = params.pot;
        this->parent->potentialReady = true;
    }

    PRISM::Multislice(params);
    {
        QMutexLocker gatekeeper(&this->parent->stackLock);
        this->parent->stack = params.stack;
        this->parent->stackReady = true;
    }

}

PotentialThread::~PotentialThread(){}
SMatrixThread::~SMatrixThread(){}
FullPRISMCalcThread::~FullPRISMCalcThread(){}
FullMultisliceCalcThread::~FullMultisliceCalcThread(){}