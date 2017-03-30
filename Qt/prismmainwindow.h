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
    friend class FullPRISMCalcThread;
    friend class FullMultisliceCalcThread;

public:
    explicit PRISMMainWindow(QWidget *parent = 0);
    ~PRISMMainWindow();

public slots:
	void setInterpolationFactor();
	void setFilenameAtoms_fromDialog();
	void setFilenameOutput_fromLineEdit();
	void setFilenameOutput_fromDialog();
//	void launch();
    void setNumGPUs(const int& numGPUs);
    void setNumThreads(const int& numThreads);
    void setNumFP(const int& numFP);
    void setPixelSize_fromLineEdit();
    void setPotBound_fromLineEdit();
    void setAlphaBeamMax_fromLineEdit();
    void setSliceThickness_fromLineEdit();
    void setCellDimX_fromLineEdit();
    void setCellDimY_fromLineEdit();
    void setCellDimZ_fromLineEdit();
    void setE0_fromLineEdit();
	void setAlgo_PRISM();
	void setAlgo_Multislice();
    void calculatePotential();
    void calculateAll();
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
    void updateSlider_lineEdits_min(int);
    void updateSlider_lineEdits_max(int);
	void updateSlider_lineEdits_max_ang(int val);
	void updateSlider_lineEdits_min_ang(int val);
    void resizeEvent(QResizeEvent* event);
    void redrawImages();
    void saveCurrentOutputImage();


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
//    QGraphicsScene *potentialScene;
    PRISM::Metadata<PRISM_FLOAT_PRECISION>* getMetadata(){return this->meta;}
    PRISM::Array3D<PRISM_FLOAT_PRECISION> potential;
    PRISM::Array3D<std::complex< PRISM_FLOAT_PRECISION> > Scompact;
    PRISM::Array4D<PRISM_FLOAT_PRECISION> output;
	PRISM::Array1D<PRISM_FLOAT_PRECISION> detectorAngles;
    QMutex potentialLock;
	QMutex sMatrixLock;
    QMutex outputLock;
    QMutex calculationLock;

    bool potentialReady;
    bool ScompactReady;
    bool outputReady;

    QImage potentialImage;
    QImage probeImage;
    QImage probeImage_pr;
    QImage probeImage_pk;
    QImage probeImage_mr;
    QImage probeImage_mk;
    QImage outputImage;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> potentialImage_float;
    PRISM::Array2D<PRISM_FLOAT_PRECISION> outputImage_float;

    PRISM_FLOAT_PRECISION contrast_potentialMin;
    PRISM_FLOAT_PRECISION contrast_potentialMax;
    PRISM_FLOAT_PRECISION contrast_outputMin;
    PRISM_FLOAT_PRECISION contrast_outputMax;
//    QImage potenetialImage;

};

unsigned char getUcharFromFloat(PRISM_FLOAT_PRECISION val,
                                PRISM_FLOAT_PRECISION contrast_low,
                                PRISM_FLOAT_PRECISION contrast_high);

#endif // PRISMMAINWINDOW_H
