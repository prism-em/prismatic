#include "prism_qthreads.h"
#include "prism_progressbar.h"
#include "PRISM01.h"
#include <iostream>

PotentialThread::PotentialThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
   //this->meta.interpolationFactor = 23;
   // std::cout << "parent interp = " << (*parent->getMetadata()).interpolationFactor << std::endl;
   // std::cout << "my interp = " << this->meta.interpolationFactor << std::endl;

}


void PotentialThread::run(){
    std::cout << "Potential thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta);
    //prism_progressbar *progress = new prism_progressbar(this->parent);
    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
    //std::cout<<"before copy this->parent->pot.at(0,0,0) = " << this->parent->pot.at(0,0,0) << std::endl;
    QMutexLocker gatekeeper(&this->parent->potentialLock);
    this->parent->potential = params.pot;
    this->parent->potentialReady = true;
    std::cout<<"after copy this->parent->pot.at(0,0,0) = " << this->parent->potential.at(0,0,0) << std::endl;

//    emit potentialReady(params.pot);
}

PotentialThread::~PotentialThread(){}


FullCalcThread::FullCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
}


void FullCalcThread::run(){
    //NEED TO FIX BY MAKING PARAMS ACCESSIBLE
    std::cout << "Full Calculation thread running" << std::endl;
    //PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta);
    PRISM::configure(meta);
    PRISM::execute_plan(meta);
//    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
//    {
//        QMutexLocker gatekeeper(&this->parent->potentialLock);
//        this->parent->potential = params.pot;
//        this->parent->potentialReady = true;
//    }

//    {
//        QMutexLocker gatekeeper(&this->parent->stackLock);
//        this->parent->stack = params.stack;
//        this->parent->stackReady = true;
//    }

//    std::cout<<"after copy this->parent->pot.at(0,0,0) = " << this->parent->potential.at(0,0,0) << std::endl;
}

FullCalcThread::~FullCalcThread(){}
