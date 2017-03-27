#include "prism_qthreads.h"
#include "PRISM01.h"
#include <iostream>

PotentialThread::PotentialThread(PRISMMainWindow *_parent) :
parent(_parent){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    this->meta = *(parent->getMetadata());
   //this->meta.interpolationFactor = 23;
   // std::cout << "parent interp = " << (*parent->getMetadata()).interpolationFactor << std::endl;
   // std::cout << "my interp = " << this->meta.interpolationFactor << std::endl;

}


void PotentialThread::run(){
    std::cout << "Potential thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta);
    PRISM::PRISM01(params);
    std::cout <<"Potential Calculated" << std::endl;
    //std::cout<<"before copy this->parent->pot.at(0,0,0) = " << this->parent->pot.at(0,0,0) << std::endl;
    this->parent->potentialLock.lock();
    this->parent->pot = params.pot;
    this->parent->potentialReady = true;
    this->parent->potentialLock.unlock();
    std::cout<<"after copy this->parent->pot.at(0,0,0) = " << this->parent->pot.at(0,0,0) << std::endl;

//    emit potentialReady(params.pot);
}
