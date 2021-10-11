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

#ifndef PRISMMAINWINDOW_H
#define PRISMMAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QMutex>
#include <QMutexLocker>
#include <QImage>
#include <QPixmap>
#include <iostream>
#include <complex>
#include "meta.h"
#include "defines.h"
#include "ArrayND.h"
#include "params.h"
#include "prism_colormapper.h"
namespace Ui
{
class PRISMMainWindow;
}

// forward declare the thread classes that run the work
class PotentialThread;

class PRISMMainWindow : public QMainWindow
{
    Q_OBJECT

    // declare the thread classes as friends so that they can
    // access the protected mutex locks and arrays
    friend class PRISMThread;
    friend class PotentialThread;
    //    friend class SMatrixThread;
    friend class ProbeThread;
    friend class FullCalcThread;
    friend class prism_progressbar;

public:
    explicit PRISMMainWindow(QWidget *parent = 0);
    bool potentialIsReady();
    bool overwriteFile();
    bool SMatrixIsReady();
    bool OutputIsReady();
    bool checkoutputArrayExists();
    bool checkoutputArrayExists_HRTEM();
    bool checkpotentialArrayExists();
    void updateUCdims(const std::string &filename);
    ~PRISMMainWindow();

public slots:
    void updateDisplay();
    void selectParameterFile();
    void writeParameterFile();
    void readParams(const std::string param_filename);
    void readAberrationFile();
    void readProbeFile();
    void setInterpolationFactorX();
    void setInterpolationFactorY();
    void setFilenameAtoms_fromDialog();
    void setFilenameOutput_fromLineEdit();
    void setFilenameOutput_fromDialog();
    void setNumGPUs(const int &numGPUs);
    void setNumThreads(const int &numThreads);
    void setNumFP(const int &numFP);
    void setNumNS(const int &numSlices);
    void setzSampling(const int &num);
    void setNumStreams(const int &numFP);
    void setPixelSizeX_fromLineEdit();
    void setPixelSizeY_fromLineEdit();
    void setPotBound_fromLineEdit();
    void setprobeSemiangle_fromLineEdit();
    void setzStart_fromLineEdit();
    // void setalphaBeamMax_fromLineEdit();
    void set2D_innerAngle_fromLineEdit();
    void set2D_outerAngle_fromLineEdit();
    void setSliceThickness_fromLineEdit();
    void setCellDimX_fromLineEdit();
    void setCellDimY_fromLineEdit();
    void setCellDimZ_fromLineEdit();
    void setTileX_fromLineEdit();
    void setTileY_fromLineEdit();
    void setTileZ_fromLineEdit();
    void setBatchGPU_fromLineEdit();
    void setBatchCPU_fromLineEdit();
    void setE0_fromLineEdit();
    void setprobeStepX_fromLineEdit();
    void setprobeStepY_fromLineEdit();
    void setAlgo_PRISM();
    void setAlgo_Multislice();
    void calculatePotential();
    void calculateAll();
    void calculateAllHRTEM();
    void calculateProbe();
    void updatePotentialImage();
    void updatePotentialDisplay();
    void updatePotentialFloatImage();
    void updateOutputImage();
    void updateOutputDisplay();
    void updateOutputFloatImage();
    void updateOutputImage_HRTEM();
    void updateOutputDisplay_HRTEM();
    void updateOutputFloatImage_HRTEM();
    void updateSliders_fromLineEdits();
    void updateSliders_fromLineEdits_ang();
    void updateContrastPotMin();
    void updateContrastPotMax();
    void updateContrastAngMin();
    void updateContrastAngMax();
    void updateSlider_lineEdits_min(int);
    void updateSlider_lineEdits_max(int);
    void updateSlider_lineEdits_max_ang(int val);
    void updateSlider_lineEdits_min_ang(int val);
    void updateAlphaMax();
    void resizeEvent(QResizeEvent *event);
    void redrawImages();
    void saveCurrentOutputImage();
    void setStreamingMode(int);
    //void toggleSaveProjectedPotential();
    void enableOutputWidgets();
    void setprobe_defocus_fromLineEdit();
    void set_dfr_min_fromLineEdit();
    void set_dfr_max_fromLineEdit();
    void set_dfr_step_fromLineEdit();
    void set_dfs_fromLineEdit();
    void setRandomSeed_fromLineEdit();
    void setprobe_C3_fromLineEdit();
    void setprobe_C5_fromLineEdit();
    void setdetectorAngleStep_fromLineEdit();
    void setprobe_Xtilt_fromLineEdit();
    void setprobe_Ytilt_fromLineEdit();
    void setxtt_min_fromLineEdit();
    void setxtt_max_fromLineEdit();
    void setxtt_step_fromLineEdit();
    void setytt_min_fromLineEdit();
    void setytt_max_fromLineEdit();
    void setytt_step_fromLineEdit();
    void setrtt_min_fromLineEdit();
    void setrtt_max_fromLineEdit();
    void setxtilt_offset_fromLineEdit();
    void setytilt_offset_fromLineEdit();
    void toggle2DOutput();
    void toggle3DOutput();
    void toggle4DOutput();
    void toggle4Dcrop();
    void set4Damax_fromLineEdit();
    void toggleDPC_CoM();
    void togglePotentialSlices();
    void togglematrixRefocus();
    void togglePotential3D();
    void toggleThermalEffects();
    void toggleSMatrixoutput();
    void toggleComplexoutput();
    void toggleProbeOutput();
    //void toggleOccupancy();
    void toggleNyquist();
    void setscan_WindowXMin_fromLineEdit();
    void setscan_WindowXMax_fromLineEdit();
    void setscan_WindowYMin_fromLineEdit();
    void setscan_WindowYMax_fromLineEdit();
    void resetCalculation();
    void resetPotential();
    void newRandomSeed();

    //    void updateProbeK_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    //    void updateProbeR_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    //    void updateProbeK_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    //    void updateProbeR_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void probeK_PRISMReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void probeR_PRISMReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void probeK_MultisliceReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void probeR_MultisliceReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void updateProbeK_PRISMDisplay();
    void updateProbeR_PRISMDisplay();
    void updateProbeK_MultisliceDisplay();
    void updateProbeR_MultisliceDisplay();
    void updateProbe_diffRDisplay();
    void updateProbe_diffKDisplay();
    void updateProbeImages();
    void calculateProbe_diffR(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void calculateProbe_diffK(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void updateAllImages();
    void update_pearsonReal(QString str);
    void update_pearsonK(QString str);
    void update_RReal(QString str);
    void update_RK(QString str);
    void potentialReceived(Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>);
    void outputReceived(Prismatic::Array4D<PRISMATIC_FLOAT_PRECISION>);
    void outputReceived_HRTEM(Prismatic::Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>>);
    void displayErrorReadingParamsDialog();
    void displayErrorReadingAtomsDialog();
    void setscan_WindowYMin_edited();
    void setscan_WindowYMax_edited();
    void setinterpYSet_edited();
    void setpixelSizeYSet_edited();
    void setprobeStepYSet_edited();
    void setprobeTiltYSet_edited();
    void checkInput_lineEdit_scanWindowXMin();
    void checkInput_lineEdit_scanWindowXMax();
    void checkInput_lineEdit_scanWindowYMin();
    void checkInput_lineEdit_scanWindowYMax();
    void checkInput_lineEdit_cellDimX();
    void checkInput_lineEdit_cellDimY();
    void checkInput_lineEdit_cellDimZ();
    void checkInput_lineEdit_tileX();
    void checkInput_lineEdit_tileY();
    void checkInput_lineEdit_tileZ();
    void checkInput_lineEdit_pixelSizeX();
    void checkInput_lineEdit_pixelSizeY();
    void checkInput_lineEdit_interpFactor_x();
    void checkInput_lineEdit_interpFactor_y();
    void userHasSetCellDims();
    void resetLinks();
    void moveBothPotentialSliders(int);
    void moveBothDetectorSliders(int);
    void openSaveAtomsDialog();
    void saveAtomCoords(QString, QString);
    void changeColormap(QString);
    bool checkProbesCalculated();
    void preventOverwrite();
    void flipOverwrite();

    //Collapse functions for each box
    void collapseSample();
    void collapseSimulation();
    void collapseStem();
    void collapseHrtem();
    void collapseOutput();
    void collapseComputational();


    //Change themes functions
    void lightField();
    void darkField();

protected:
    void setFilenameAtoms(const std::string &filename);
    void setFilenameOutput(const std::string &filename);
    void setRealspacePixelSize(const PRISMATIC_FLOAT_PRECISION &pixel_size);
    void setPotBound(const PRISMATIC_FLOAT_PRECISION &potBound);
    void setNumFP(const size_t &numFP);
    void setNumNS(const size_t &numSlices);
    void setE0(const PRISMATIC_FLOAT_PRECISION &E0);
    // void setAlphaBeamMax(const PRISMATIC_FLOAT_PRECISION &alphaBeamMax);
    void setSliceThickness(const PRISMATIC_FLOAT_PRECISION &thickness);
    void setCellDimX(const int &dimX);
    void setCellDimY(const int &dimY);
    void setCellDimZ(const int &dimZ);
    void setAlgo(const Prismatic::Algorithm algo);

private:
    Ui::PRISMMainWindow *ui;
    Prismatic::Colormapper colormapper;
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> *meta;
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars;
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars_multi;
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> *getMetadata() { return this->meta; }
    Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION> potential;
    Prismatic::Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> smatrix;
    Prismatic::Array4D<PRISMATIC_FLOAT_PRECISION> output;
    Prismatic::Array1D<PRISMATIC_FLOAT_PRECISION> detectorAngles;
    std::vector<PRISMATIC_FLOAT_PRECISION> pixelSize;

    QMutex potentialLock;
    QMutex outputLock;
    QMutex dataLock;
    QMutex calculationLock;
    QMutex probeLock;

    bool potentialReady;
    bool ScompactReady;
    bool outputReady;
    //bool saveProjectedPotential;
    bool probeSetupReady;
    bool potentialArrayExists;
    bool outputArrayExists;
    bool outputArrayExists_HRTEM;
    bool probesCalculated;
    bool interpYSet;
    bool pixelSizeYSet;
    bool probeStepYSet;
    bool probeTiltYSet;
    bool minWindowYSet;
    bool maxWindowYSet;
    bool overwriteCheck = false;


    //collapse funtction variables
    bool sampleClosed = true;
    bool simulationClosed = true;
    bool stemClosed = true;
    bool hrtemClosed = true;
    bool outputClosed = true;
    bool computationalClosed = true;

    //Height each box gets
    //Computational and Output boxes both have unique sizes
    int boxOpen = 260;
    int boxClosed = 20;
    int scrollOpen = 230;

    int animSpeed = 600;
    
    QImage potentialImage;
    QImage probeImage;
    QImage probeImage_pr;
    QImage probeImage_pk;
    QImage probeImage_mr;
    QImage probeImage_mk;
    QImage probeImage_diffr;
    QImage probeImage_diffk;
    QImage outputImage;
    QImage outputImage_HRTEM;

    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> potentialImage_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> outputImage_float; 
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> outputImage_HRTEM_float; 
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> probeImage_pr_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> probeImage_pk_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> probeImage_mr_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> probeImage_mk_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> probeImage_diffr_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> probeImage_diffk_float;

    PRISMATIC_FLOAT_PRECISION contrast_potentialMin;
    PRISMATIC_FLOAT_PRECISION contrast_potentialMax;
    PRISMATIC_FLOAT_PRECISION contrast_outputMin;
    PRISMATIC_FLOAT_PRECISION contrast_outputMax;
    PRISMATIC_FLOAT_PRECISION contrast_outputMin_HRTEM;
    PRISMATIC_FLOAT_PRECISION contrast_outputMax_HRTEM;
    PRISMATIC_FLOAT_PRECISION currently_calculated_X;
    PRISMATIC_FLOAT_PRECISION currently_calculated_Y;
};

unsigned char getUcharFromFloat(PRISMATIC_FLOAT_PRECISION val,
                                PRISMATIC_FLOAT_PRECISION contrast_low,
                                PRISMATIC_FLOAT_PRECISION contrast_high);

#endif // PRISMMAINWINDOW_H
