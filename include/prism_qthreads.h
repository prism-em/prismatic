#ifndef PRISM_QTHREADS_H
#define PRISM_QTHREADS_H
#include <QThread>
#include "ArrayND.h"
#include "params.h"
#include "defines.h"
#include "prismmainwindow.h"
// defines QThread derived classes for running work from the PRISM GUI
//class PRISMMainWindow;

class PotentialThread : public QThread {
    Q_OBJECT
    void run() Q_DECL_OVERRIDE;
    friend class PRISMMainWindow;
public:
    PotentialThread(PRISMMainWindow *_parent);


//public slots:
   // signals:
   // void potentialReady(PRISM::Array3D<PRISM_FLOAT_PRECISION>);

private:
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
};
#endif // PRISM_QTHREADS_H
