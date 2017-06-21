// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_QTHREADS_H
#define PRISM_QTHREADS_H
#include <QThread>
#include "prism_progressbar.h"
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
    explicit PotentialThread(PRISMMainWindow *_parent, prism_progressbar *progressbar);
    virtual ~PotentialThread();
//signals:
//    void potentialCalculated();
private:
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
};


class SMatrixThread : public QThread {
    Q_OBJECT
    void run() Q_DECL_OVERRIDE;
    friend class PRISMMainWindow;
public:
    explicit SMatrixThread(PRISMMainWindow *_parent, prism_progressbar *progressbar);
    virtual ~SMatrixThread();
//signals:
//    void ScompactCalculated();
private:
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
};

class FullCalcThread : public QThread {
    Q_OBJECT
    void run() Q_DECL_OVERRIDE;
    friend class PRISMMainWindow;
public:
    explicit FullCalcThread(PRISMMainWindow *_parent, prism_progressbar *progressbar);
    virtual ~FullCalcThread();
signals:
    void potentialCalculated();
    void ScompactCalculated();
    void stackCalculated();
private:
    Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
};

#endif // PRISM_QTHREADS_H
