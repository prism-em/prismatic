// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

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
    friend class PotentialThread;
    friend class SMatrixThread;
    friend class ProbeThread;
    friend class FullPRISMCalcThread;
    friend class FullMultisliceCalcThread;
    friend class prism_progressbar;

public:
    explicit PRISMMainWindow(QWidget *parent = 0);
	bool potentialIsReady();
	bool SMatrixIsReady();
	bool OutputIsReady();
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
    void setPixelSize_fromLineEdit();
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
    void setprobe_stepX_fromLineEdit();
    void setprobe_stepY_fromLineEdit();
	void setAlgo_PRISM();
	void setAlgo_Multislice();
    void calculatePotential();
    void calculateSMatrix();
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
    void toggleStreamingMode();
    void toggleSaveProjectedPotential();
    void enableOutputWidgets();
    void setprobe_defocus_fromLineEdit();
    void setRandomSeed_fromLineEdit();
    void setprobe_C3_fromLineEdit();
    void setprobe_C5_fromLineEdit();
    void setDetector_angle_step_fromLineEdit();
    void setprobe_Xtilt_fromLineEdit();
    void setprobe_Ytilt_fromLineEdit();
    void toggle3DOutput();
    void toggle4DOutput();
    void toggleThermalEffects();
    void setscan_WindowXMin_fromLineEdit();
    void setscan_WindowXMax_fromLineEdit();
    void setscan_WindowYMin_fromLineEdit();
    void setscan_WindowYMax_fromLineEdit();
    void resetCalculation();
    void newRandomSeed();
    void updateProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void updateProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void updateProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void updateProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void updateProbe_diffR(PRISM::Array2D<PRISM_FLOAT_PRECISION>, PRISM::Array2D<PRISM_FLOAT_PRECISION> arr_contrast);
    void updateProbe_diffK(PRISM::Array2D<PRISM_FLOAT_PRECISION>, PRISM::Array2D<PRISM_FLOAT_PRECISION> arr_contrast);
    void update_pearsonReal(QString str);
    void update_pearsonK(QString str);
    void update_RReal(QString str);
    void update_RK(QString str);
//signals:
//    void signalProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
//    void signalProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
//    void signalProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
//    void signalProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
//    void testImage();

protected:

    void setFilenameAtoms(const std::string& filename);
    void setFilenameOutput(const std::string& filename);
    void setRealspacePixelSize(const PRISM_FLOAT_PRECISION& pixel_size);
    void setPotBound(const PRISM_FLOAT_PRECISION& potBound);
    void setNumFP(const size_t& numFP);
    void setE0(const PRISM_FLOAT_PRECISION& E0);
    void setAlphaBeamMax(const PRISM_FLOAT_PRECISION& alphaBeamMax);
    void setSliceThickness(const PRISM_FLOAT_PRECISION& thickness);
    void setCellDimX(const int& dimX);
    void setCellDimY(const int& dimY);
    void setCellDimZ(const int& dimZ);
	void setAlgo(const PRISM::Algorithm algo);

private:
    Ui::PRISMMainWindow *ui;
    PRISM::Metadata<PRISM_FLOAT_PRECISION> *meta;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> pars;
    PRISM::Parameters<PRISM_FLOAT_PRECISION> pars_multi;

//    QGraphicsScene *potentialScene;
    PRISM::Metadata<PRISM_FLOAT_PRECISION>* getMetadata(){return this->meta;}
//    PRISM::Array3D<PRISM_FLOAT_PRECISION> potential;
//    PRISM::Array3D<std::complex< PRISM_FLOAT_PRECISION> > Scompact;
//    PRISM::Array3D<PRISM_FLOAT_PRECISION> output;
    PRISM::Array1D<PRISM_FLOAT_PRECISION> detectorAngles;

//    QMutex potentialLock;
//    QMutex sMatrixLock;
//    QMutex outputLock;
    QMutex dataLock;
    QMutex calculationLock;

    bool potentialReady;
    bool ScompactReady;
    bool outputReady;
    bool saveProjectedPotential;
    bool probeSetupReady;
    bool potentialImageExists;
    bool outputImageExists;

    QImage potentialImage;
    QImage probeImage;
    QImage probeImage_pr;
    QImage probeImage_pk;
    QImage probeImage_mr;
    QImage probeImage_mk;
    QImage probeImage_diffr;
    QImage probeImage_diffk;
    QImage outputImage;

    PRISM::Array2D<PRISM_FLOAT_PRECISION> potentialImage_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> outputImage_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> probeImage_pr_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> probeImage_pk_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> probeImage_mr_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> probeImage_mk_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> probeImage_diffr_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> probeImage_diffk_float;

    PRISM_FLOAT_PRECISION contrast_potentialMin;
    PRISM_FLOAT_PRECISION contrast_potentialMax;
    PRISM_FLOAT_PRECISION contrast_outputMin;
    PRISM_FLOAT_PRECISION contrast_outputMax;
    PRISM_FLOAT_PRECISION currently_calculated_X;
    PRISM_FLOAT_PRECISION currently_calculated_Y;
//    QImage potenetialImage;

};

unsigned char getUcharFromFloat(PRISM_FLOAT_PRECISION val,
                                PRISM_FLOAT_PRECISION contrast_low,
                                PRISM_FLOAT_PRECISION contrast_high);

#endif // PRISMMAINWINDOW_H
