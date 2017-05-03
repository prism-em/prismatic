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
private:
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
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
signals:
    void potentialCalculated();
    void ScompactCalculated();
private:
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
};

class ProbeThread : public QThread {
    Q_OBJECT
    void run() Q_DECL_OVERRIDE;
    friend class PRISMMainWindow;
public:
    explicit ProbeThread(PRISMMainWindow *_parent, PRISM_FLOAT_PRECISION _X, PRISM_FLOAT_PRECISION _Y, prism_progressbar *progressbar);
    virtual ~ProbeThread();
signals:
    void potentialCalculated();
    void ScompactCalculated();
    void signalProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void signalProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void signalProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void signalProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void signalProbe_diffK(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
    void signalProbe_diffR(PRISM::Array2D<PRISM_FLOAT_PRECISION>);
private:
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    PRISM_FLOAT_PRECISION X, Y;
    prism_progressbar *progressbar;
};

class FullPRISMCalcThread : public QThread {
    Q_OBJECT
    void run() Q_DECL_OVERRIDE;
    friend class PRISMMainWindow;
public:
    explicit FullPRISMCalcThread(PRISMMainWindow *_parent, prism_progressbar *progressbar);
    virtual ~FullPRISMCalcThread();
signals:
    void potentialCalculated();
	void ScompactCalculated();
	void outputCalculated();
private:
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
};

class FullMultisliceCalcThread : public QThread {
    Q_OBJECT
    void run() Q_DECL_OVERRIDE;
    friend class PRISMMainWindow;
public:
    explicit FullMultisliceCalcThread(PRISMMainWindow *_parent, prism_progressbar *progressbar);
    virtual ~FullMultisliceCalcThread();
signals:
	void potentialCalculated();
	void ScompactCalculated();
	void outputCalculated();
private:
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
};

#endif // PRISM_QTHREADS_H
