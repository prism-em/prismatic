#ifndef PRISM_PROGRESSBAR_H
#define PRISM_PROGRESSBAR_H

#include <QDialog>
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
	void setText(QString str);
public slots:
    void setStepPotential();
    void update_calculatingPotential(long, long);
//    void update_calculatingSMatrix(int, int);
//    void update_calculatingPRISM(int, int, int, int);
//    void update_calculatingMultislice(int, int, int, int);

private:
    Ui::prism_progressbar *ui;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;

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
