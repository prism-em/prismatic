#include "prism_progressbar.h"
#include "ui_prism_progressbar.h"
#include <algorithm>
prism_progressbar::prism_progressbar(PRISMMainWindow *_parent) :
    parent(_parent),
    ui(new Ui::prism_progressbar),
    potentialCurrentSlice(0),
    potentialTotalSlices(0),
    SMatrixCurrentBeam(0),
    SMatrixTotalBeams(0),
    PRISMCurrentXProbe(0),
    PRISMTotalXProbes(0),
    PRISMCurrentYProbe(0),
    PRISMTotalYProbes(0),
    MultisliceCurrentXProbe(0),
    MultisliceTotalXProbes(0),
    MultisliceCurrentYProbe(0),
    MultisliceTotalYProbes(0)
{
    ui->setupUi(this);
    ui->progressBar->setValue(0);
	connect(this, SIGNAL(updateCalcStatus(QString)), this, SLOT(updateCalcStatusMessage(QString)));
	connect(this, SIGNAL(updateProgressBar(int)), ui->progressBar, SLOT(setValue(int)));

}

void prism_progressbar::setStepPotential(){
    ui->lbl_calcStatus->setText(QString("Computing Projected Potential Slices"));
}

//void prism_progressbar::setAlgorithmPRISM(){
//    ui->lbl_calcStatus->setText("Computing Projected Potential Slices");
//}

void prism_progressbar::update_calculatingPotential(long current, long total){
    potentialCurrentSlice = std::max(potentialCurrentSlice, current);
    ui->lbl_Description->setText(QString("Slice ") +
                                 QString::number(current) +
                                 QString("/") +
                                 QString::number(total));
}

//void prism_progressbar::setText(const QString str){
//    ui->lbl_Description->setText(str);
//}
void prism_progressbar::updateCalcStatusMessage(const QString str){
    ui->lbl_calcStatus->setText(str);
}
void prism_progressbar::signalCalcStatusMessage(const QString str){
    emit updateCalcStatus(str);
}
void prism_progressbar::signalPotentialUpdate(const long current, const long total){
	std::lock_guard<std::mutex> gatekeeper(dataLock);
	potentialCurrentSlice = std::max(potentialCurrentSlice, current);
	emit updateCalcStatus(QString("Slice ") +
	                      QString::number(potentialCurrentSlice) +
	                      QString("/") +
	                      QString::number(total));
	emit updateProgressBar(100*total/current);
}
prism_progressbar::~prism_progressbar()
{
    delete ui;
}
