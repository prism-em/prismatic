#ifndef PRISM_PROGRESSBAR_H
#define PRISM_PROGRESSBAR_H

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
public slots:
    void setStepPotential();
    void update_calculatingPotential(long, long);
    void updateCalcStatusMessage(const QString str);

//    void update_calculatingSMatrix(int, int);
//    void update_calculatingPRISM(int, int, int, int);
//    void update_calculatingMultislice(int, int, int, int);
signals:
    void updateCalcStatus(const QString str);
	void updateProgressBar(int value);
private:
    Ui::prism_progressbar *ui;
    PRISMMainWindow *parent;
    //prism_progressbar *progressbar;
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
    std::mutex dataLock;
    long potentialCurrentSlice;
    long potentialTotalSlices;
    long SMatrixCurrentBeam;
    long SMatrixTotalBeams;
    long PRISMCurrentXProbe;
    long PRISMTotalXProbes;
    long PRISMCurrentYProbe;
    long PRISMTotalYProbes;
    long MultisliceCurrentXProbe;
    long MultisliceTotalXProbes;
    long MultisliceCurrentYProbe;
    long MultisliceTotalYProbes;

};

#endif // PRISM_PROGRESSBAR_H
