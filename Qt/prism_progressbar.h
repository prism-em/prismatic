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

#ifndef PRISM_PROGRESSBAR_H
#define PRISM_PROGRESSBAR_H

//#define PRISMATIC_BUILDING_GUI

#include <QDialog>
#include <mutex>
#include "prismmainwindow.h"
#include "meta.h"
#include "configure.h"
namespace Ui {
class prism_progressbar;
}

class prism_progressbar : public QDialog
{
    Q_OBJECT

public:
    explicit prism_progressbar(PRISMMainWindow *parent);
    ~prism_progressbar();
    //void setText(const QString str);
    void signalCalcStatusMessage(const QString str);
    void signalPotentialUpdate(const long current, const long total);
	void setProgress(int);
    void signalDescriptionMessage(const QString str);
    void signalScompactUpdate(const long current, const long total);
    void signalOutputUpdate(const long current, const long total);
    void resetOutputs();
public slots:
    void setStepPotential();
    void update_calculatingPotential(long, long);
    void updateCalcStatusMessage(const QString str);
    void updateDescription(const QString str);
    void setTitle(const QString str);

//    void update_calculatingSMatrix(int, int);
//    void update_calculatingPRISM(int, int, int, int);
//    void update_calculatingMultislice(int, int, int, int);
signals:
    void updateDescriptionMessage(const QString str);
    void updateCalcStatus(const QString str);
	void updateProgressBar(int value);

private:
    Ui::prism_progressbar *ui;
    PRISMMainWindow *parent;
    //prism_progressbar *progressbar;
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    std::mutex dataLock;
    long potentialCurrentSlice;
    long potentialTotalSlices;
    long SMatrixCurrentBeam;
    long SMatrixTotalBeams;
    long currentProbe;
    long totalProbes;

};

#endif // PRISM_PROGRESSBAR_H
