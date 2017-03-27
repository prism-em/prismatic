#ifndef PRISMMAINWINDOW_H
#define PRISMMAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QMutex>
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
public:
    explicit PRISMMainWindow(QWidget *parent = 0);
    ~PRISMMainWindow();
   // PRISM::Metadata<PRISM_FLOAT_PRECISION>* getMetadata(){return this->meta;}
  //  PRISM::Array3D<PRISM_FLOAT_PRECISION> pot;
  //  PRISM::Array3D<std::complex< PRISM_FLOAT_PRECISION> > Scompact;
  //  PRISM::Array4D<PRISM_FLOAT_PRECISION> stack;
  //  QMutex potentialLock;
   // QMutex stackLock;
   // QMutex ScompactLock;

public slots:
	void setInterpolationFactor();
	void setFilenameAtoms_fromDialog();
	void setFilenameOutput_fromLineEdit();
	void setFilenameOutput_fromDialog();
	void launch();
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
//    void testImage();

protected:
    Ui::PRISMMainWindow *ui;
    PRISM::Metadata<PRISM_FLOAT_PRECISION> *meta;
    QGraphicsScene *potentialScene;

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

    PRISM::Metadata<PRISM_FLOAT_PRECISION>* getMetadata(){return this->meta;}
    PRISM::Array3D<PRISM_FLOAT_PRECISION> pot;
    PRISM::Array3D<std::complex< PRISM_FLOAT_PRECISION> > Scompact;
    PRISM::Array4D<PRISM_FLOAT_PRECISION> stack;
    QMutex potentialLock;

    bool potentialReady;
    //bool ScompactReady;
    //bool stackReady;


};

#endif // PRISMMAINWINDOW_H
