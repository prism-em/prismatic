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
namespace Ui {
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
    friend class FullPRISMCalcThread;
    friend class FullMultisliceCalcThread;
    friend class prism_progressbar;

public:
    explicit PRISMMainWindow(QWidget *parent = 0);
	bool potentialIsReady();
	bool SMatrixIsReady();
	bool OutputIsReady();
    bool checkoutputArrayExists();
    bool checkpotentialArrayExists();
    void updateUCdims(const std::string& filename);
    ~PRISMMainWindow();

public slots:
	void setInterpolationFactorX();
	void setInterpolationFactorY();
	void setFilenameAtoms_fromDialog();
	void setFilenameOutput_fromLineEdit();
	void setFilenameOutput_fromDialog();
//	void launch();
    void setNumGPUs(const int& numGPUs);
    void setNumThreads(const int& numThreads);
    void setNumFP(const int& numFP);
    void setNumStreams(const int& numFP);
    void setPixelSizeX_fromLineEdit();
    void setPixelSizeY_fromLineEdit();
    void setPotBound_fromLineEdit();
    void setprobeSemiangle_fromLineEdit();
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
//    void calculateSMatrix();
    void calculateAll();
    void calculateProbe();
    void updatePotentialImage();
    void updatePotentialDisplay();
    void updatePotentialFloatImage();
    void updateOutputImage();
    void updateOutputDisplay();
    void updateOutputFloatImage();
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
    void resizeEvent(QResizeEvent* event);
    void redrawImages();
    void saveCurrentOutputImage();
    void setStreamingMode(int);
    void toggleSaveProjectedPotential();
    void enableOutputWidgets();
    void setprobe_defocus_fromLineEdit();
    void setRandomSeed_fromLineEdit();
    void setprobe_C3_fromLineEdit();
    void setprobe_C5_fromLineEdit();
    void setdetectorAngleStep_fromLineEdit();
    void setprobe_Xtilt_fromLineEdit();
    void setprobe_Ytilt_fromLineEdit();
    void toggle3DOutput();
    void toggle4DOutput();
    void toggleThermalEffects();
    void toggleOccupancy();
    void setscan_WindowXMin_fromLineEdit();
    void setscan_WindowXMax_fromLineEdit();
    void setscan_WindowYMin_fromLineEdit();
    void setscan_WindowYMax_fromLineEdit();
    void resetCalculation();
    void newRandomSeed();
    void updateProbeK_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void updateProbeR_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void updateProbeK_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void updateProbeR_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>);
    void updateProbe_diffR(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>, Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr_contrast);
    void updateProbe_diffK(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>, Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr_contrast);
    void update_pearsonReal(QString str);
    void update_pearsonK(QString str);
    void update_RReal(QString str);
    void update_RK(QString str);
    void potentialReceived(Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>);
    void outputReceived(Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>);
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
    void updateSlider_PotentialCombo(int);
    void openSaveAtomsDialog();
    void saveAtomCoords(QString, QString);
protected:
    void setFilenameAtoms(const std::string& filename);
    void setFilenameOutput(const std::string& filename);
    void setRealspacePixelSize(const PRISMATIC_FLOAT_PRECISION& pixel_size);
    void setPotBound(const PRISMATIC_FLOAT_PRECISION& potBound);
    void setNumFP(const size_t& numFP);
    void setE0(const PRISMATIC_FLOAT_PRECISION& E0);
    void setAlphaBeamMax(const PRISMATIC_FLOAT_PRECISION& alphaBeamMax);
    void setSliceThickness(const PRISMATIC_FLOAT_PRECISION& thickness);
    void setCellDimX(const int& dimX);
    void setCellDimY(const int& dimY);
    void setCellDimZ(const int& dimZ);
	void setAlgo(const Prismatic::Algorithm algo);

private:
    Ui::PRISMMainWindow *ui;
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> *meta;
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars;
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars_multi;

//    QGraphicsScene *potentialScene;
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION>* getMetadata(){return this->meta;}
    Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION> potential;
//    Prismatic::Array3D<std::complex< PRISMATIC_FLOAT_PRECISION> > Scompact;
    Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION> output;
    Prismatic::Array1D<PRISMATIC_FLOAT_PRECISION> detectorAngles;
    std::vector<PRISMATIC_FLOAT_PRECISION> pixelSize;

    QMutex potentialLock;
//    QMutex sMatrixLock;
    QMutex outputLock;
    QMutex dataLock;
    QMutex calculationLock;

    bool potentialReady;
    bool ScompactReady;
    bool outputReady;
    bool saveProjectedPotential;
    bool probeSetupReady;
    bool potentialArrayExists;
    bool outputArrayExists;
    bool interpYSet;
    bool pixelSizeYSet;
    bool probeStepYSet;
    bool probeTiltYSet;
    bool minWindowYSet;
    bool maxWindowYSet;

    QImage potentialImage;
    QImage probeImage;
    QImage probeImage_pr;
    QImage probeImage_pk;
    QImage probeImage_mr;
    QImage probeImage_mk;
    QImage probeImage_diffr;
    QImage probeImage_diffk;
    QImage outputImage;

    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> potentialImage_float;
    Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> outputImage_float;
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
    PRISMATIC_FLOAT_PRECISION currently_calculated_X;
    PRISMATIC_FLOAT_PRECISION currently_calculated_Y;
//    QImage potenetialImage;

};

unsigned char getUcharFromFloat(PRISMATIC_FLOAT_PRECISION val,
                                PRISMATIC_FLOAT_PRECISION contrast_low,
                                PRISMATIC_FLOAT_PRECISION contrast_high);

#endif // PRISMMAINWINDOW_H
