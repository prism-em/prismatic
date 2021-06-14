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

#include "prism_qthreads.h"
#include "prism_progressbar.h"
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include "PRISM03_calcOutput.h"
#include "Multislice_calcOutput.h"
#include "Multislice_entry.h"
#include "PRISM_entry.h"
#include "HRTEM_entry.h"
#include "utility.h"
#include <iostream>
#include <string>
#include "H5Cpp.h"
#include "QMessageBox"
#include <stdio.h>
#include "fileIO.h"
#include "aberration.h"

PRISMThread::PRISMThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) : parent(_parent), progressbar(_progressbar)
{
    // construct the thread with a copy of the metadata so that any upstream changes don't mess with this calculation
    QMutexLocker gatekeeper(&this->parent->dataLock);
    this->meta = *(parent->getMetadata());
}

void PRISMThread::passPotentialToParent(Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION> &arr)
{
    this->parent->potentialReceived(arr);
    emit potentialCalculated();
}

PotentialThread::PotentialThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) : PRISMThread(_parent, _progressbar){};

void PotentialThread::run()
{
    // create parameters
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> params(meta, this, progressbar);
    //prevent it from trying to write to a non-existing H5 file
    params.meta.savePotentialSlices = false;
    progressbar->signalDescriptionMessage("Computing projected potential");

    {
        QMutexLocker calculationLocker(&this->parent->calculationLock);
        Prismatic::configure(meta);
        // calculate potential
        Prismatic::PRISM01_calcPotential(params);
        this->parent->potentialReceived(params.pot);
        emit potentialCalculated();
    }
    // acquire the mutex so we can safely copy to the GUI copy of the potential
    QMutexLocker gatekeeper(&this->parent->dataLock);
    // perform copy
    this->parent->pars = params;

    // indicate that the potential is ready
    std::cout << "Projected potential calculation complete" << std::endl;
}

ProbeThread::ProbeThread(PRISMMainWindow *_parent, PRISMATIC_FLOAT_PRECISION _X, PRISMATIC_FLOAT_PRECISION _Y, prism_progressbar *_progressbar, bool _use_log_scale) : PRISMThread(_parent, _progressbar), X(_X), Y(_Y), use_log_scale(_use_log_scale){};

void ProbeThread::run()
{
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> params(meta, this, progressbar);
    progressbar->signalDescriptionMessage("Computing single probe");
    progressbar->setProgress(10);
    //    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> params(meta);

    //    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> params_multi(params);
    params.meta.save4DOutput = false;
    params.meta.savePotentialSlices = false;
    params.meta.saveDPC_CoM = false;

    QMutexLocker calculationLocker(&this->parent->calculationLock);

    if ((!this->parent->potentialIsReady()) || !(params.meta == *(this->parent->getMetadata())))
    {
        //    if ((!this->parent->potentialIsReady())){

        Prismatic::configure(meta);
        this->parent->resetCalculation(); // any time we are computing the potential we are effectively starting over the whole calculation, so make sure all flags are reset
        Prismatic::PRISM01_calcPotential(params);
        std::cout << "Potential Calculated" << std::endl;
        {
            QMutexLocker gatekeeper(&this->parent->dataLock);
            this->parent->pars = params;
        }
        //            this->parent->potentialArrayExists = true;

        this->parent->potentialReceived(params.pot);
        emit potentialCalculated();
        /*
        if (this->parent->saveProjectedPotential){
          std::string outfile = Prismatic::remove_extension(params.meta.filenameOutput.c_str());
          outfile += ("_potential.mrc");
          std::cout<<"Outputting potential with filename "<<outfile<<std::endl;
          params.pot.toMRC_f(outfile.c_str());
        }
        */
    }
    else
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "Potential already calculated. Using existing result." << std::endl;
    }

    if ((!this->parent->SMatrixIsReady()) || !(params.meta == *(this->parent->getMetadata())))
    {
        Prismatic::PRISM02_calcSMatrix(params);
        {
            std::cout << "S-Matrix finished calculating." << std::endl;
            QMutexLocker gatekeeper(&this->parent->dataLock);
            // perform copy
            this->parent->pars = params;

            // indicate that the S-Matrix is ready
            this->parent->ScompactReady = true;
        }
    }
    else
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        std::cout << "S-Matrix already calculated. Using existing result." << std::endl;
    }

    std::pair<Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>>, Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>>> prism_probes, multislice_probes;
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> params_multi(params);
    params_multi.meta.save4DOutput = false;
    params_multi.meta.savePotentialSlices = false;
    params_multi.meta.saveDPC_CoM = false;

    // setup and calculate PRISM probe
    Prismatic::setupCoordinates_2(params);

    // setup angles of detector and image sizes
    Prismatic::setupDetector(params);

    // setup coordinates and indices for the beams
    Prismatic::setupBeams_2(params);

    // setup Fourier coordinates for the S-matrix
    Prismatic::setupFourierCoordinates(params);

    // initialize the output to the correct size for the output mode
    Prismatic::createStack_integrate(params);

    //  perform some necessary setup transformations of the data
    Prismatic::transformIndices(params);

    // initialize/compute the probes
    Prismatic::initializeProbes(params);

    std::cout << "Getting PRISM Probe" << std::endl;
    std::cout << "X = " << X << std::endl;
    std::cout << "Y = " << Y << std::endl;
    prism_probes = Prismatic::getSinglePRISMProbe_CPU(params, X, Y);

    // setup and calculate Multislice probe
    // setup coordinates and build propagators
    Prismatic::setupCoordinates_multislice(params_multi);

    // setup detector coordinates and angles
    Prismatic::setupDetector_multislice(params_multi);

    // create initial probes
    Prismatic::setupProbes_multislice(params_multi);

    // create transmission array
    Prismatic::createTransmission(params_multi);

    // initialize output stack
    Prismatic::createStack(params_multi);

    multislice_probes = Prismatic::getSingleMultisliceProbe_CPU(params_multi, X, Y);

    //    if (params.meta == *(this->parent->getMetadata())){
    //        std::cout << "Metadata equal!" << std::endl;
    //    } else {
    //        std::cout << "Metadata not equal!" << std::endl;
    //    }

    QMutexLocker gatekeeper(&this->parent->dataLock);
    // perform copy
    //    this->parent->pars = params;
    this->parent->probeSetupReady = true;

    prism_probes = upsamplePRISMProbe(prism_probes.first,
                                      multislice_probes.first.get_dimj(),
                                      multislice_probes.first.get_dimi(),
                                      Y / params.pixelSize[0] / 2,
                                      X / params.pixelSize[1] / 2);

    PRISMATIC_FLOAT_PRECISION pr_sum, pk_sum, mr_sum, mk_sum;
    pr_sum = pk_sum = mr_sum = mk_sum = 0;
    for (auto &i : prism_probes.first)
        pr_sum += abs(i);
    for (auto &i : prism_probes.second)
        pk_sum += abs(i);
    for (auto &i : multislice_probes.first)
        mr_sum += abs(i);
    for (auto &i : multislice_probes.second)
        mk_sum += abs(i);

    for (auto &i : prism_probes.first)
        i /= pr_sum;
    for (auto &i : prism_probes.second)
        i /= pk_sum;
    for (auto &i : multislice_probes.first)
        i /= mr_sum;
    for (auto &i : multislice_probes.second)
        i /= mk_sum;

    emit signal_pearsonReal(QString("Pearson Correlation = ") + QString::number(computePearsonCorrelation(prism_probes.first, multislice_probes.first)));
    emit signal_pearsonK(QString("Pearson Correlation = ") + QString::number(computePearsonCorrelation(prism_probes.second, multislice_probes.second)));
    emit signal_RReal(QString("R = ") + QString::number(computeRfactor(prism_probes.first, multislice_probes.first)));
    emit signal_RK(QString("R = ") + QString::number(computeRfactor(prism_probes.second, multislice_probes.second)));

    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> pr = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{prism_probes.first.get_dimj(), prism_probes.first.get_dimi()}});
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> pk = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{prism_probes.second.get_dimj(), prism_probes.second.get_dimi()}});
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> mr = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> mk = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{multislice_probes.second.get_dimj(), multislice_probes.second.get_dimi()}});
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> diffr = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{multislice_probes.first.get_dimj(), multislice_probes.first.get_dimi()}});
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> diffk = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{multislice_probes.second.get_dimj(), multislice_probes.second.get_dimi()}});

    for (auto i = 0; i < prism_probes.first.size(); ++i)
    {
        pr[i] = std::abs(prism_probes.first[i]);
    }
    for (auto i = 0; i < prism_probes.second.size(); ++i)
    {
        pk[i] = std::abs(prism_probes.second[i]);
    }
    for (auto i = 0; i < multislice_probes.first.size(); ++i)
    {
        mr[i] = std::abs(multislice_probes.first[i]);
    }
    for (auto i = 0; i < multislice_probes.second.size(); ++i)
    {
        mk[i] = std::abs(multislice_probes.second[i]);
    }
    for (auto i = 0; i < prism_probes.second.size(); ++i)
    {
        diffr[i] = (std::abs(pr[i] - mr[i]));
        diffk[i] = (std::abs(pk[i] - mk[i]));
    }

    emit signalProbeR_PRISM((pr));
    emit signalProbeK_PRISM(fftshift2(pk));
    emit signalProbeR_Multislice((mr));
    emit signalProbeK_Multislice(fftshift2(mk));
    emit signalProbe_diffR(diffr);
    emit signalProbe_diffK(fftshift2(diffk));
}

FullCalcThread::FullCalcThread(PRISMMainWindow *_parent, prism_progressbar *_progressbar) : PRISMThread(_parent, _progressbar){};

void FullCalcThread::run()
{

    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> params(meta, this, progressbar);
    if (!Prismatic::testFilenameOutput(params.meta.filenameOutput.c_str()))
    {
        std::cout << "Aborting calculation, please choose an accessible output directory" << std::endl;
        return;
    }

    if (Prismatic::testFilenameOutput(params.meta.filenameOutput.c_str()) == 2)
    {
        emit overwriteWarning();
        if (this->parent->overwriteFile())
        {
            //params.outputFile.flush(H5F_SCOPE_GLOBAL);
            remove(params.meta.filenameOutput.c_str());
            this->thread()->sleep(1);
            this->parent->flipOverwrite(); //flip the check back
        }
        else
        {
            return;
        }
    }

    progressbar->signalDescriptionMessage("Initiating Multislice simulation");

    std::cout << "Also do CPU work: " << params.meta.alsoDoCPUWork << std::endl;
    QMutexLocker calculationLocker(&this->parent->calculationLock);
    Prismatic::configure(meta);

    if( (!this->parent->potentialIsReady()) || !(meta == *(this->parent->getMetadata())) ){
        this->parent->resetCalculation();
    }
    else
    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        params = this->parent->pars;
        params.progressbar = progressbar;
        params.parent_thread = this;
        std::cout << "Potential already calculated. Using existing result for first frozen phonon. " << std::endl;
    }
    
    params.potentialReady = this->parent->potentialIsReady();
    if(params.meta.algorithm == Prismatic::Algorithm::Multislice){
        Prismatic::Multislice_entry_pars(params);
    }
    else if (params.meta.algorithm == Prismatic::Algorithm::PRISM){
        Prismatic::PRISM_entry_pars(params);
    }
    else if (params.meta.algorithm == Prismatic::Algorithm::HRTEM)
    {
        Prismatic::HRTEM_entry_pars(params);
    }

    {
        QMutexLocker gatekeeper(&this->parent->dataLock);
        this->parent->pars = params;
        gatekeeper.unlock();
    }

    {
        QMutexLocker gatekeeper(&this->parent->outputLock);
        this->parent->detectorAngles = params.detectorAngles;
        for (auto &a : this->parent->detectorAngles)
            a *= 1000; // convert to mrads
        this->parent->pixelSize = params.pixelSize;

        gatekeeper.unlock();
    }


    if(params.meta.algorithm == Prismatic::Algorithm::HRTEM){
        this->parent->outputReceived_HRTEM(params.Scompact);
    }
    else{
        this->parent->outputReceived(params.output);
    }
    emit outputCalculated();
    std::cout << "Calculation complete" << std::endl;

    calculationLocker.unlock();
}

PRISMThread::~PRISMThread() {}
PotentialThread::~PotentialThread() {}
ProbeThread::~ProbeThread() {}
FullCalcThread::~FullCalcThread() {}
