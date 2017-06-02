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
    QMutexLocker gatekeeper(&this->parent->dataLock);
    this->meta = *(parent->getMetadata());
}

void PotentialThread::run(){
    // create parameters
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
    {
        QMutexLocker calculationLocker(&this->parent->calculationLock);
        PRISM::configure(meta);
        // calculate potential
        PRISM::PRISM01(params);
        this->parent->potentialReceived(params.pot);
        emit potentialCalculated();
    }
    // acquire the mutex so we can safely copy to the GUI copy of the potential
    QMutexLocker gatekeeper(&this->parent->dataLock);
    // perform copy
    this->parent->pars = params;
    // indicate that the potential is ready

//    this->parent->potentialArrayExists = true;
    if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
    std::cout << "Projected potential calculation complete" << std::endl;
}

//SMatrixThread::SMatrixThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
//        parent(_parent), progressbar(_progressbar){
//    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
//    QMutexLocker gatekeeper(&this->parent->dataLock);
//    this->meta = *(parent->getMetadata());
//}

//void SMatrixThread::run(){
//    QMutexLocker calculationLocker(&this->parent->calculationLock);
//    PRISM::configure(meta);
//    // create parameters
//    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
//    // calculate potential if it hasn't been already
//    if (!this->parent->potentialIsReady()){
//        // calculate potential
//        PRISM::PRISM01(params);
//	    QMutexLocker gatekeeper(&this->parent->dataLock);
//        // indicate that the potential is ready

////        this->parent->potentialArrayExists = true;
//        if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
//        this->parent->potentialReceived(params.pot);
//        emit potentialCalculated();
//    } else {
//        QMutexLocker gatekeeper(&this->parent->dataLock);
//        params = this->parent->pars;
//        params.progressbar = progressbar;
//        std::cout << "Potential already calculated. Using existing result." << std::endl;
//    }


////    // calculate S-Matrix
//    PRISM::PRISM02(params);
//    QMutexLocker gatekeeper(&this->parent->dataLock);

////    // perform copy
//    this->parent->pars = params;
//    this->parent->ScompactReady = true;
//    std::cout << "S-matrix calculation complete" << std::endl;
//}


ProbeThread::ProbeThread(PRISMMainWindow *_parent, PRISM_FLOAT_PRECISION _X, PRISM_FLOAT_PRECISION _Y, prism_progressbar *_progressbar, bool _use_log_scale) :
parent(_parent), X(_X), Y(_Y), progressbar(_progressbar), use_log_scale(_use_log_scale){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    QMutexLocker gatekeeper(&this->parent->dataLock);
    this->meta = *(parent->getMetadata());
}

void ProbeThread::run(){
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
//    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta);

//    PRISM::Parameters<PRISM_FLOAT_PRECISION> params_multi(params);

    QMutexLocker calculationLocker(&this->parent->calculationLock);

    if (!this->parent->potentialIsReady()){
        PRISM::configure(meta);
        this->parent->resetCalculation(); // any time we are computing the potential we are effectively starting over the whole calculation, so make sure all flags are reset
        PRISM::PRISM01(params);
        std::cout <<"Potential Calculated" << std::endl;
        {
            QMutexLocker gatekeeper(&this->parent->dataLock);
            this->parent->pars = params;
        }
//            this->parent->potentialArrayExists = true;
        if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
        this->parent->potentialReceived(params.pot);
        emit potentialCalculated();
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "Potential already calculated. Using existing result." << std::endl;
    }

    if (!this->parent->SMatrixIsReady()){
        PRISM::PRISM02(params);
    {
//        std::cout << "S-Matrix finished calculating." << std::endl;
//        QMutexLocker gatekeeper(&this->parent->dataLock);
//        // perform copy
//        this->parent->pars = params;

//        // indicate that the S-Matrix is ready
//        this->parent->ScompactReady = true;
    }
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
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

    std::cout << "Getting PRISM Probe" << std::endl;
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
    std::cout << "Getting Multislice Probe" << std::endl;

    multislice_probes = PRISM::getSingleMultisliceProbe_CPU(params_multi, X, Y);


    QMutexLocker gatekeeper(&this->parent->dataLock);
    // perform copy
//    this->parent->pars = params;
    this->parent->probeSetupReady = true;

    prism_probes = upsamplePRISMProbe(prism_probes.first,
                                      multislice_probes.first.get_dimj(),
                                      multislice_probes.first.get_dimi(),
                                      Y / params.pixelSize[0] / 2 ,
                                      X / params.pixelSize[1] / 2);

    PRISM_FLOAT_PRECISION pr_sum, pk_sum, mr_sum, mk_sum;
    pr_sum = pk_sum  = mr_sum = mk_sum = 0;
    for (auto& i : prism_probes.first)       pr_sum += abs(i);
    for (auto& i : prism_probes.second)      pk_sum += abs(i);
    for (auto& i : multislice_probes.first)  mr_sum += abs(i);
    for (auto& i : multislice_probes.second) mk_sum += abs(i);

    for (auto& i : prism_probes.first)       i /= pr_sum;
    for (auto& i : prism_probes.second)      i /= pk_sum;
    for (auto& i : multislice_probes.first)  i /= mr_sum;
    for (auto& i : multislice_probes.second) i /= mk_sum;

    std::cout << "emitting signal" << std::endl;
    emit signal_pearsonReal(QString("Pearson Correlation = ") + QString::number(computePearsonCorrelation(prism_probes.first, multislice_probes.first)));
    emit signal_pearsonK(QString("Pearson Correlation = ") + QString::number(computePearsonCorrelation(prism_probes.second, multislice_probes.second)));
    emit signal_RReal(QString("R = ") + QString::number(computeRfactor(prism_probes.first, multislice_probes.first)));
    emit signal_RK(QString("R = ") + QString::number(computeRfactor(prism_probes.second, multislice_probes.second)));
//    pr_sum = std::accumulate(&prism_probes.first[0], &*(prism_probes.first.end()-1), (PRISM_FLOAT_PRECISION)0.0, [](std::complex<PRISM_FLOAT_PRECISION> a){return abs(a);});

//    PRISM::Array2D<PRISM_FLOAT_PRECISION> debug = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
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

//    for (auto j = 0; j < prism_probes.first.get_dimj(); ++j){
//        for (auto i = 0; i < prism_probes.first.get_dimi(); ++i){
//        debug.at(j,i) = std::abs(prism_probes.first.at(j,i));
//        }
//    }
//    debug.toMRC_f("/mnt/spareA/clion/PRISM/build/db_p.mrc");
//    for (auto j = 0; j < prism_probes.second.get_dimj(); ++j){
//        for (auto i = 0; i < prism_probes.second.get_dimi(); ++i){
//        debug.at(j,i) = std::abs(prism_probes.second.at(j,i));
//        }
//    }
//    debug.toMRC_f("/mnt/spareA/clion/PRISM/build/dbk_p.mrc");


//    std::cout << "prism_probes.first.get_dimj() = " << prism_probes.first.get_dimj() <<std::endl;
//    std::cout << "multislice_probes.first.get_dimj() = " << multislice_probes.first.get_dimj() <<std::endl;

    PRISM::Array2D<PRISM_FLOAT_PRECISION> pr    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{prism_probes.first.get_dimj(), prism_probes.first.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> pk    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{prism_probes.second.get_dimj(), prism_probes.second.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> mr    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> mk    = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.second.get_dimj(), multislice_probes.second.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> diffr = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    PRISM::Array2D<PRISM_FLOAT_PRECISION> diffk = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{multislice_probes.second.get_dimj(), multislice_probes.second.get_dimi()}});



if (use_log_scale){
    for (auto i = 0; i < prism_probes.first.size(); ++i){
        pr[i] =  std::log(1e-5 + std::abs(prism_probes.first[i]));
    }
    for (auto i = 0; i < prism_probes.second.size(); ++i){
        pk[i] =  std::log(1e-5 + std::abs(prism_probes.second[i]));
    }
    for (auto i = 0; i < multislice_probes.first.size(); ++i){
        mr[i] =  std::log(1e-5 + std::abs(multislice_probes.first[i]));
    }
    for (auto i = 0; i < multislice_probes.second.size(); ++i){
        mk[i] =  std::log(1e-5 + std::abs(multislice_probes.second[i]));
    }
    for (auto i = 0; i < prism_probes.second.size(); ++i){
        diffr[i] =  log(1e-5 + std::abs(pr[i] - mr[i]));
        diffk[i] =  log(1e-5 + std::abs(pk[i] - mk[i]));
    }
} else{
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
    for (auto i = 0; i < prism_probes.second.size(); ++i){
        diffr[i] =  (std::abs(pr[i] - mr[i]));
        diffk[i] =  (std::abs(pk[i] - mk[i]));
    }
}

    emit signalProbeR_PRISM((pr));
    emit signalProbeK_PRISM(fftshift2(pk));
    emit signalProbeR_Multislice((mr));
    emit signalProbeK_Multislice(fftshift2(mk));
    emit signalProbe_diffR((diffr), mr);
    emit signalProbe_diffK(fftshift2(diffk), mk);
}

FullPRISMCalcThread::FullPRISMCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    QMutexLocker gatekeeper(&this->parent->dataLock);
    this->meta = *(parent->getMetadata());
}


void FullPRISMCalcThread::run(){
    std::cout << "Full PRISM Calculation thread running" << std::endl;
    emit signalTitle("PRISM: Frozen Phonon #1");
    bool error_reading = false;
    QMutexLocker gatekeeper(&this->parent->dataLock);
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params;
    try {
        params = PRISM::Parameters<PRISM_FLOAT_PRECISION>(meta,progressbar);
    }catch (...){
        std::cout <<"An error occurred while attempting to read from file " << meta.filename_atoms << std::endl;
        error_reading = true;
    }
    gatekeeper.unlock();


    if (error_reading){
        emit signalErrorReadingAtomsDialog();
    } else {
    QMutexLocker calculationLocker(&this->parent->calculationLock);

    PRISM::configure(meta);
//  //  PRISM::Parameters<PRISM_FLOAT_PRECISION> params = PRISM::execute_plan(meta);
    if (!this->parent->potentialIsReady()){
        this->parent->resetCalculation(); // any time we are computing the potential we are effectively starting over the whole calculation, so make sure all flags are reset
        PRISM::PRISM01(params);
        std::cout <<"Potential Calculated" << std::endl;
    {
//        QMutexLocker gatekeeper(&this->parent->potentialLock);
        QMutexLocker gatekeeper(&this->parent->dataLock);
        this->parent->pars = params;
//        this->parent->potential = params.pot;

//        this->parent->potentialArrayExists = true;
        if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
    }
    this->parent->potentialReceived(params.pot);
    emit potentialCalculated();
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "Potential already calculated. Using existing result." << std::endl;
    }

    if (!this->parent->SMatrixIsReady()){
    PRISM::PRISM02(params);
    {
//        QMutexLocker gatekeeper(&this->parent->dataLock);

//        // perform copy
//        this->parent->pars = params;

//        // indicate that the potential is ready
//        this->parent->ScompactReady = true;
    }
    } else {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "S-Matrix already calculated. Using existing result." << std::endl;
    }
//    emit ScompactCalculated();

    PRISM::PRISM03(params);



    if (params.meta.numFP > 1) {
        // run the rest of the frozen phonons
        PRISM::Array3D<PRISM_FLOAT_PRECISION> net_output(params.output);
        for (auto fp_num = 1; fp_num < params.meta.numFP; ++fp_num){
            PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
            params.meta.random_seed = rand() % 100000;
            emit signalTitle("PRISM: Frozen Phonon #" + QString::number(1 + fp_num));
            progressbar->resetOutputs();
            PRISM::PRISM01(params);
            PRISM::PRISM02(params);
            PRISM::PRISM03(params);
            net_output += params.output;
        }
        // divide to take average
        for (auto&i:net_output) i/=params.meta.numFP;

        params.output = net_output;
    }



//    {
//        QMutexLocker gatekeeper(&this->parent->dataLock);
//        this->parent->pars = params;
//    }

    {
        QMutexLocker gatekeeper(&this->parent->outputLock);
        this->parent->detectorAngles = params.detectorAngles;
        for (auto& a:this->parent->detectorAngles) a*=1000; // convert to mrads
        this->parent->pixelSize = params.pixelSize;

//        this->parent->outputArrayExists = true;

//        params.output.toMRC_f(params.meta.filename_output.c_str());
    }
	if (params.meta.save3DOutput)params.output.toMRC_f(params.meta.filename_output.c_str());
//    this->parent->outputReceived(params.output);
    this->parent->outputReceived(params.output);
    emit outputCalculated();
    std::cout << "PRISM calculation complete" << std::endl;
    }
}



FullMultisliceCalcThread::FullMultisliceCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) :
        parent(_parent), progressbar(_progressbar){
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    QMutexLocker gatekeeper(&this->parent->dataLock);
    this->meta = *(parent->getMetadata());
}


void FullMultisliceCalcThread::run(){
    std::cout << "Full Multislice Calculation thread running" << std::endl;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
	QMutexLocker calculationLocker(&this->parent->calculationLock);
    PRISM::configure(meta);
        if (!this->parent->potentialIsReady()) {
            this->parent->resetCalculation(); // any time we are computing the potential we are effectively starting over the whole calculation, so make sure all flags are reset
            PRISM::PRISM01(params);
            std::cout << "Potential Calculated" << std::endl;
            {
                QMutexLocker gatekeeper(&this->parent->dataLock);
                this->parent->pars = params;
                //        this->parent->potential = params.pot;

//                this->parent->potentialArrayExists = true;
                if (this->parent->saveProjectedPotential)params.pot.toMRC_f("potential.mrc");
            }
        } else {
            QMutexLocker gatekeeper(&this->parent->dataLock);
            params = this->parent->pars;
            params.progressbar = progressbar;
            std::cout << "Potential already calculated. Using existing result." << std::endl;
        }

    this->parent->potentialReceived(params.pot);
    emit potentialCalculated();

    PRISM::Multislice(params);


    if (params.meta.numFP > 1) {
        // run the rest of the frozen phonons
        PRISM::Array3D<PRISM_FLOAT_PRECISION> net_output(params.output);
        for (auto fp_num = 1; fp_num < params.meta.numFP; ++fp_num){
            PRISM::Parameters<PRISM_FLOAT_PRECISION> params(meta, progressbar);
            params.meta.random_seed = rand() % 100000;
            emit signalTitle("PRISM: Frozen Phonon #" + QString::number(1 + fp_num));
            progressbar->resetOutputs();
            PRISM::PRISM01(params);
            PRISM::Multislice(params);
            net_output += params.output;
        }
        // divide to take average
        for (auto&i:net_output) i/=params.meta.numFP;
        params.output = net_output;
    }


    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        this->parent->pars = params;
    }
    {
        QMutexLocker gatekeeper(&this->parent->outputLock);
	    this->parent->detectorAngles = params.detectorAngles;
	    for (auto& a:this->parent->detectorAngles) a*=1000; // convert to mrads
        this->parent->pixelSize = params.pixelSize;

//        this->parent->outputArrayExists = true;
//        PRISM::Array3D<PRISM_FLOAT_PRECISION> reshaped_output = PRISM::zeros_ND<3, PRISM_FLOAT_PRECISION>(
//        {{params.output.get_diml(), params.output.get_dimk(), params.output.get_dimj()}});
//        auto ptr = reshaped_output.begin();
//        for (auto &i:params.output)*ptr++=i;
//        reshaped_output.toMRC_f(params.meta.filename_output.c_str());
//        params.output.toMRC_f(params.meta.filename_output.c_str());
    }
	if (params.meta.save3DOutput)params.output.toMRC_f(params.meta.filename_output.c_str());
    this->parent->outputReceived(params.output);
    emit outputCalculated();
    std::cout << "Multislice calculation complete" << std::endl;
}

PotentialThread::~PotentialThread(){}
//SMatrixThread::~SMatrixThread(){}
ProbeThread::~ProbeThread(){}
FullPRISMCalcThread::~FullPRISMCalcThread(){}
FullMultisliceCalcThread::~FullMultisliceCalcThread(){}
