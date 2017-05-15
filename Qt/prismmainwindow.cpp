// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include <QFileDialog>
#include "prismmainwindow.h"
#include "ui_prismmainwindow.h"
#include <fstream>
#include <iostream>
#include <utility>
//#include "PRISM_entry.h"
#include "configure.h"
#include "prism_qthreads.h"
#include "prism_progressbar.h"
#include "utility.h"

bool validateFilename(const std::string str){
    std::ifstream f(str);
    return f.good();
}
PRISM_FLOAT_PRECISION calculateLambda(PRISM::Metadata<PRISM_FLOAT_PRECISION> meta);

PRISMMainWindow::PRISMMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PRISMMainWindow),
    potentialReady(false),
    ScompactReady(false),
    outputReady(false),
    saveProjectedPotential(false),
    probeSetupReady(false),
    potentialImageExists(false),
    outputImageExists(false),
    potentialImage(QImage()),
    currently_calculated_X(0.0),
    currently_calculated_Y(0.0)
{
    qRegisterMetaType<PRISM::Array2D< PRISM_FLOAT_PRECISION> >("PRISM::Array2D<PRISM_FLOAT_PRECISION>");
	// build Qt generated interface
    ui->setupUi(this);

	// set window title
    setWindowTitle("PRISM");

    ui->box_sampleSettings->setStyleSheet("QGroupBox { \
                                          border: 1px solid gray;\
                                          border-radius: 9px;\
                                          margin-top: 0.5em;\
                                      }  QGroupBox::title {\
                                         subcontrol-origin: margin;\
                                         left: 150px;\
                                         padding: 0 3px 0 3px;\
                                         }");

//ui->box_calculationSettings->setStyleSheet("QGroupBox { \
//                                      border: 1px solid gray;\
//                                      border-radius: 9px;\
//                                      margin-top: 0.5em;\
//                                  }  QGroupBox::title {\
//                                     font-weight: bold;\
//                                     subcontrol-origin: margin;\
//                                     left: 145px;\
//                                     padding: 0 300px 0 3px;\
//                                     }");

ui->box_calculationSettings->setStyleSheet("QGroupBox { \
                                      border: 1px solid gray;\
                                      border-radius: 9px;\
                                      margin-top: 0.5em;\
                                  }  QGroupBox::title {\
                                     font-weight: bold;\
                                     subcontrol-origin: margin;\
                                     left: 145px;\
                                     padding: 0 3px 0 3px;\
                                     }");


    ui->box_simulationSettings->setStyleSheet("QGroupBox { \
                                        border: 1px solid gray;\
                                        border-radius: 9px;\
                                        margin-top: 0.5em;\
                                    }  QGroupBox::title {\
                                       subcontrol-origin: margin;\
                                       left: 150px;\
                                       padding: 0 3px 0px 3px;\
                                       }");


    potentialImage.load(":/images/prism.png");
    probeImage.load(":/images/probe.png");
    outputImage.load(":/images/output.png");
    probeImage_pr.load(":/images/airy.png");
    probeImage_pk.load(":/images/airy.png");
    probeImage_mr.load(":/images/airy.png");
    probeImage_mk.load(":/images/airy.png");
    probeImage_diffk.load(":/images/airy.png");
    probeImage_diffr.load(":/images/airy.png");
    redrawImages();
    ui->lbl_image_potential->setPixmap(QPixmap::fromImage(potentialImage.scaled(ui->lbl_image_potential->width(),
                                                                                ui->lbl_image_potential->height(),
                                                                                Qt::KeepAspectRatio)));
	// set initially displayed values based on the default parameters
	this->meta = new PRISM::Metadata<PRISM_FLOAT_PRECISION>;
	{
		std::stringstream ss;
        ss << this->meta->interpolationFactorX;
        this->ui->lineedit_interpFactor_x->setText(QString::fromStdString(ss.str()));
		ss.str("");
        ss << this->meta->interpolationFactorY;
        this->ui->lineedit_interpFactor_y->setText(QString::fromStdString(ss.str()));
        ss.str("");
		ss << this->meta->potBound;
		this->ui->lineEdit_potbound->setText(QString::fromStdString(ss.str()));
		ss.str("");
        ss << (this->meta->probeSemiangle * 1e3);
        this->ui->lineEdit_probeSemiangle->setText(QString::fromStdString(ss.str()));
		ss.str("");
		ss << this->meta->sliceThickness;
		this->ui->lineEdit_sliceThickness->setText(QString::fromStdString(ss.str()));
		ss.str("");
		ss << this->meta->cellDim[2];
		this->ui->lineEdit_cellDimX->setText(QString::fromStdString(ss.str()));
		ss.str("");
		ss << this->meta->cellDim[1];
		this->ui->lineEdit_cellDimY->setText(QString::fromStdString(ss.str()));
		ss.str("");
		ss << this->meta->cellDim[0];
		this->ui->lineEdit_cellDimZ->setText(QString::fromStdString(ss.str()));
		ss.str("");
        ss << this->meta->tileX;
        this->ui->lineEdit_tileX->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << this->meta->tileY;
        this->ui->lineEdit_tileY->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << this->meta->tileZ;
        this->ui->lineEdit_tileZ->setText(QString::fromStdString(ss.str()));
        ss.str("");
		ss << (this->meta->E0 *1e-3);
		this->ui->lineEdit_E0->setText(QString::fromStdString(ss.str()));
		ss.str("");
        ss << (this->meta->probe_stepX);
        this->ui->lineEdit_probeStepX->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->probe_stepY);
        this->ui->lineEdit_probeStepY->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->probeDefocus * 1e10);
        this->ui->lineEdit_probeDefocus->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->probeXtilt);
        this->ui->lineEdit_probeTiltX->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->probeYtilt);
        this->ui->lineEdit_probeTiltY->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->detector_angle_step * 1e3);
        this->ui->lineEdit_detectorAngle->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->scanWindowXMin);
        this->ui->lineEdit_scanWindowXMin->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->scanWindowXMax);
        this->ui->lineEdit_scanWindowXMax->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->scanWindowYMin);
        this->ui->lineEdit_scanWindowYMin->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->scanWindowYMax);
        this->ui->lineEdit_scanWindowYMax->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->random_seed);
        this->ui->lineEdit_randomSeed->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->batch_size_CPU);
        this->ui->lineEdit_batchCPU->setText(QString::fromStdString(ss.str()));
        ss.str("");
        ss << (this->meta->batch_size_GPU);
        this->ui->lineEdit_batchGPU->setText(QString::fromStdString(ss.str()));
        ss.str("");


		this->ui->lineedit_outputfile->setText(QString::fromStdString(ss.str()));
		this->ui->spinBox_numGPUs->setValue(this->meta->NUM_GPUS);
		this->ui->spinBox_numThreads->setValue(this->meta->NUM_THREADS);
		this->ui->spinBox_numFP->setValue(this->meta->numFP);

	}

    ui->checkBox_thermalEffects->setChecked(meta->include_thermal_effects);
    ui->checkBox_3D->setChecked(meta->save3DOutput);
    ui->checkBox_4D->setChecked(meta->save4DOutput);

	switch (this->meta->algorithm){
		case PRISM::Algorithm::PRISM :      this->ui->radBtn_PRISM->setChecked(true);
			                                this->ui->radBtn_Multislice->setChecked(false);
			break;
		case PRISM::Algorithm::Multislice : this->ui->radBtn_PRISM->setChecked(false);
										    this->ui->radBtn_Multislice->setChecked(true);
			break;
	}
#ifndef PRISM_ENABLE_GPU
	this->ui->spinBox_numGPUs->setEnabled(false);
    this->ui->checkBox_streamdata->setEnabled(false);
#endif //PRISM_ENABLE_GPU


    ui->lbl_angstrom->setText(QString::fromUtf8("\u212B"));
    ui->lbl_sliceThickness->setText(QString::fromUtf8("Slice\nThickness (\u212B)"));
    ui->lbl_probeStep->setText(QString::fromUtf8("Probe\nStep (\u212B)"));
//    updateAlphaMax();
    ui->lbl_alphaMax->setText(QString::fromUtf8("\u03B1 max = ??"));
    ui->lbl_lambda->setText(QString::fromUtf8("\u03BB = ") + QString::number(calculateLambda(*meta)) + QString::fromUtf8("\u212B"));
    ui->lbl_potBound->setText(QString::fromUtf8("Potential\nBound (\u212B)"));
    ui->lbl_pixelSize->setText(QString::fromUtf8("Pixel\nSize (\u212B)"));
    ui->lbl_defocus->setText(QString::fromUtf8("C1 (defocus)(\u212B)"));
    ui->label_Xprobe->setText(QString::fromUtf8("X (\u212B)"));
    ui->label_Yprobe->setText(QString::fromUtf8("Y (\u212B)"));

    this->ui->lineedit_outputfile->setText(QString::fromStdString(this->meta->filename_output));

    // connect signals and slots
    connect(this->ui->lineedit_interpFactor_x,SIGNAL(editingFinished()),this,SLOT(setInterpolationFactorX()));
    connect(this->ui->lineedit_interpFactor_y,SIGNAL(editingFinished()),this,SLOT(setInterpolationFactorY()));
	connect(this->ui->lineedit_outputfile,SIGNAL(editingFinished()),this,SLOT(setFilenameOutput_fromLineEdit()));
	connect(this->ui->btn_atomsfile_browse, SIGNAL(pressed()), this, SLOT(setFilenameAtoms_fromDialog()));
    connect(this->ui->spinBox_numGPUs, SIGNAL(valueChanged(int)), this, SLOT(setNumGPUs(const int&)));
    connect(this->ui->spinBox_numThreads, SIGNAL(valueChanged(int)), this, SLOT(setNumThreads(const int&)));
    connect(this->ui->spinBox_numFP, SIGNAL(valueChanged(int)), this, SLOT(setNumFP(const int&)));
    connect(this->ui->lineEdit_probeSemiangle, SIGNAL(editingFinished()), this, SLOT(setprobeSemiangle_fromLineEdit()));
    connect(this->ui->lineedit_pixelSize, SIGNAL(editingFinished()), this, SLOT(setPixelSize_fromLineEdit()));
    connect(this->ui->lineEdit_batchCPU, SIGNAL(editingFinished()), this, SLOT(setBatchCPU_fromLineEdit()));
    connect(this->ui->lineEdit_batchGPU, SIGNAL(editingFinished()), this, SLOT(setBatchGPU_fromLineEdit()));
    connect(this->ui->lineEdit_potbound, SIGNAL(editingFinished()), this, SLOT(setPotBound_fromLineEdit()));
    connect(this->ui->lineEdit_sliceThickness, SIGNAL(editingFinished()), this, SLOT(setSliceThickness_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimX, SIGNAL(editingFinished()), this, SLOT(setCellDimX_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimY, SIGNAL(editingFinished()), this, SLOT(setCellDimY_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimZ, SIGNAL(editingFinished()), this, SLOT(setCellDimZ_fromLineEdit()));
    connect(this->ui->lineEdit_tileX, SIGNAL(editingFinished()), this, SLOT(setTileX_fromLineEdit()));
    connect(this->ui->lineEdit_tileY, SIGNAL(editingFinished()), this, SLOT(setTileY_fromLineEdit()));
    connect(this->ui->lineEdit_tileZ, SIGNAL(editingFinished()), this, SLOT(setTileZ_fromLineEdit()));
    connect(this->ui->lineEdit_randomSeed, SIGNAL(editingFinished()), this, SLOT(setRandomSeed_fromLineEdit()));
    connect(this->ui->lineEdit_probeDefocus, SIGNAL(editingFinished()), this, SLOT(setprobe_defocus_fromLineEdit()));
    connect(this->ui->lineEdit_C3, SIGNAL(editingFinished()), this, SLOT(setprobe_C3_fromLineEdit()));
    connect(this->ui->lineEdit_C5, SIGNAL(editingFinished()), this, SLOT(setprobe_C5_fromLineEdit()));
    connect(this->ui->lineEdit_detectorAngle, SIGNAL(editingFinished()), this, SLOT(setDetector_angle_step_fromLineEdit()));
    connect(this->ui->lineEdit_probeTiltX, SIGNAL(editingFinished()), this, SLOT(setprobe_Xtilt_fromLineEdit()));
    connect(this->ui->lineEdit_probeTiltY, SIGNAL(editingFinished()), this, SLOT(setprobe_Ytilt_fromLineEdit()));
    connect(this->ui->lineEdit_probeStepX, SIGNAL(editingFinished()), this, SLOT(setprobe_stepX_fromLineEdit()));
    connect(this->ui->lineEdit_probeStepY, SIGNAL(editingFinished()), this, SLOT(setprobe_stepY_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowXMin, SIGNAL(editingFinished()), this, SLOT(setscan_WindowXMin_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowXMax, SIGNAL(editingFinished()), this, SLOT(setscan_WindowXMax_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowYMin, SIGNAL(editingFinished()), this, SLOT(setscan_WindowYMin_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowYMax, SIGNAL(editingFinished()), this, SLOT(setscan_WindowYMax_fromLineEdit()));
    connect(this->ui->lineEdit_E0, SIGNAL(editingFinished()), this, SLOT(setE0_fromLineEdit()));
	connect(this->ui->radBtn_PRISM, SIGNAL(clicked(bool)), this, SLOT(setAlgo_PRISM()));
	connect(this->ui->radBtn_Multislice, SIGNAL(clicked(bool)), this, SLOT(setAlgo_Multislice()));
    connect(this->ui->btn_calcPotential, SIGNAL(clicked(bool)), this, SLOT(calculatePotential()));
//    connect(this->ui->btn_calcSmatrix, SIGNAL(clicked(bool)), this, SLOT(calculateSMatrix()));
    connect(this->ui->btn_go, SIGNAL(clicked(bool)), this, SLOT(calculateAll()));
    connect(this->ui->lineEdit_slicemin, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits()));
    connect(this->ui->lineEdit_slicemax, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits()));
    connect(this->ui->slider_slicemin, SIGNAL(valueChanged(int)), this, SLOT(updateSlider_lineEdits_min(int)));
    connect(this->ui->slider_slicemax, SIGNAL(valueChanged(int)), this, SLOT(updateSlider_lineEdits_max(int)));
    connect(this->ui->slider_slicemin, SIGNAL(valueChanged(int)), this, SLOT(updatePotentialFloatImage()));
    connect(this->ui->slider_slicemax, SIGNAL(valueChanged(int)), this, SLOT(updatePotentialFloatImage()));
	connect(this->ui->slider_angmin, SIGNAL(valueChanged(int)), this, SLOT(updateSlider_lineEdits_min_ang(int)));
	connect(this->ui->slider_angmax, SIGNAL(valueChanged(int)), this, SLOT(updateSlider_lineEdits_max_ang(int)));
	connect(this->ui->slider_angmin, SIGNAL(valueChanged(int)), this, SLOT(updateOutputFloatImage()));
	connect(this->ui->slider_angmax, SIGNAL(valueChanged(int)), this, SLOT(updateOutputFloatImage()));
	connect(this->ui->lineEdit_angmin, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits_ang()));
	connect(this->ui->lineEdit_angmax, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits_ang()));
//    connect(this->ui->lineEdit_angmin, SIGNAL(textChanged(QString)), this, SLOT(updateSliders_fromLineEdits_ang()));
//    connect(this->ui->lineEdit_angmax, SIGNAL(textChanged(QString)), this, SLOT(updateSliders_fromLineEdits_ang()));
    connect(this->ui->contrast_outputMin, SIGNAL(editingFinished()), this, SLOT(updateContrastAngMin()));
    connect(this->ui->contrast_outputMax, SIGNAL(editingFinished()), this, SLOT(updateContrastAngMax()));
    connect(this->ui->lineEdit_contrastPotMin, SIGNAL(editingFinished()), this, SLOT(updateContrastPotMin()));
    connect(this->ui->lineEdit_contrastPotMax, SIGNAL(editingFinished()), this, SLOT(updateContrastPotMax()));
    connect(this->ui->tabs, SIGNAL(currentChanged(int)), this, SLOT(redrawImages()));
    connect(this->ui->btn_saveOutputImage, SIGNAL(clicked(bool)), this, SLOT(saveCurrentOutputImage()));
    connect(this->ui->checkBox_streamdata, SIGNAL(toggled(bool)), this, SLOT(toggleStreamingMode()));
    connect(this->ui->checkBox_saveProjectedPotential, SIGNAL(toggled(bool)), this, SLOT(toggleSaveProjectedPotential()));
    connect(this->ui->btn_reset, SIGNAL(clicked()), this, SLOT(resetCalculation()));
    connect(this->ui->btn_calculateProbe, SIGNAL(clicked()), this, SLOT(calculateProbe()));
    connect(this->ui->checkBox_3D, SIGNAL(toggled(bool)), this, SLOT(toggle3DOutput()));
    connect(this->ui->checkBox_4D, SIGNAL(toggled(bool)), this, SLOT(toggle4DOutput()));
    connect(this->ui->checkBox_thermalEffects, SIGNAL(toggled(bool)), this, SLOT(toggleThermalEffects()));

    this->ui->checkBox_streamdata->setChecked(false);
     updateAlphaMax();
}



void PRISMMainWindow::setAlgo_PRISM(){
	std::cout << "Setting algorithm to PRISM" << std::endl;
	setAlgo(PRISM::Algorithm::PRISM);
}
void PRISMMainWindow::setAlgo_Multislice(){
	std::cout << "Setting algorithm to Multislice" << std::endl;
	setAlgo(PRISM::Algorithm::Multislice);
}
void PRISMMainWindow::setAlgo(const PRISM::Algorithm algo){
	this->meta->algorithm = algo;
    resetCalculation();
}

void PRISMMainWindow::setInterpolationFactorX(){
	bool flag;
    const size_t& new_f = this->ui->lineedit_interpFactor_x->text().toUInt(&flag);
	if (flag){
		std::cout << "Setting interpolation factor X to " << new_f << std::endl;
		this->meta->interpolationFactorX = new_f;
	} else{
        std::cout << "Invalid interpolation factor X input: " <<  this->ui->lineedit_interpFactor_x->text().toStdString() << std::endl;
	}
    resetCalculation();
}

void PRISMMainWindow::setInterpolationFactorY(){
    bool flag;
    const size_t& new_f = this->ui->lineedit_interpFactor_x->text().toUInt(&flag);
    if (flag){
        std::cout << "Setting interpolation factor Y to " << new_f << std::endl;
        this->meta->interpolationFactorY = new_f;
    } else{
        std::cout << "Invalid interpolation factor Y input: " <<  this->ui->lineedit_interpFactor_x->text().toStdString() << std::endl;
    }
    resetCalculation();
}


void PRISMMainWindow::setFilenameAtoms_fromDialog(){
	QString filename;
    filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename, tr("Atomic Model(*.csv *.CSV *.xyz *.XYZ);;All files(*)"));
    if (validateFilename(filename.toStdString())){
        this->setFilenameAtoms(filename.toStdString());
        ui->btn_go->setEnabled(true);
        ui->btn_calcPotential->setEnabled(true);
        this->setWindowTitle(QString::fromStdString(std::string("PRISM (") + std::string(filename.toStdString() + std::string(")"))));
    }
    resetCalculation();
}

void PRISMMainWindow::setFilenameOutput_fromDialog(){
	QString filename;
	filename = QFileDialog::getOpenFileName(this, tr("AnyFile"), filename);
	this->ui->lineedit_outputfile->setText(filename);
	this->setFilenameOutput(filename.toStdString());
    resetCalculation();
}


void PRISMMainWindow::setFilenameOutput_fromLineEdit(){
	const std::string& filename = this->ui->lineedit_outputfile->text().toStdString();
	this->setFilenameOutput(filename);
    resetCalculation();
}

void PRISMMainWindow::setFilenameAtoms(const std::string& filename){
	std::cout << "Setting atoms filename to " << filename << std::endl;
	this->meta->filename_atoms = filename;
    resetCalculation();
}

void PRISMMainWindow::setFilenameOutput(const std::string& filename){
	std::cout << "Setting output filename to " << filename << std::endl;
	this->meta->filename_output = filename;
    resetCalculation();
}


void PRISMMainWindow::setNumGPUs(const int& num){
    if (num >= 0){
        this->meta->NUM_GPUS = num;
        std::cout << "Setting number of GPUs to " << num << std::endl;
    }
}

void PRISMMainWindow::setNumThreads(const int& num){
    if (num > 0){
        this->meta->NUM_THREADS = num;
        std::cout << "Setting number of CPU Threads to " << num << std::endl;
    }
}

void PRISMMainWindow::setNumFP(const int& num){
    if (num > 0){
        this->meta->numFP = num;
        std::cout << "Setting number of frozen phonon configurations to " << num << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setPixelSize_fromLineEdit(){
    PRISM_FLOAT_PRECISION val =(PRISM_FLOAT_PRECISION)this->ui->lineedit_pixelSize->text().toDouble();
    if (val > 0){
        this->meta->realspace_pixelSize = val;
        this->meta->realspace_pixelSize = val;
        std::cout << "Setting X/Y pixel size to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setPotBound_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_potbound->text().toDouble();
    if (val > 0){
        this->meta->potBound = val;
        std::cout << "Setting potential bound to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobeSemiangle_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeSemiangle->text().toDouble();
    if (val > 0){
        this->meta->probeSemiangle = val / 1000;
        std::cout << "Setting probe semiangle to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setSliceThickness_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_sliceThickness->text().toDouble();
    if (val > 0){
        this->meta->sliceThickness = val;
        std::cout << "Setting sliceThickness to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setCellDimX_fromLineEdit(){
    double val = this->ui->lineEdit_cellDimX->text().toDouble();
    if (val > 0){
        this->meta->cellDim[2] = (PRISM_FLOAT_PRECISION)val;
        std::cout << "Setting X cell dimension to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setCellDimY_fromLineEdit(){
    double val = this->ui->lineEdit_cellDimY->text().toDouble();
    if (val > 0){
        this->meta->cellDim[1] = (PRISM_FLOAT_PRECISION)val;
        std::cout << "Setting Y cell dimension to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setCellDimZ_fromLineEdit(){
    double val = this->ui->lineEdit_cellDimZ->text().toDouble();
    if (val > 0){
        this->meta->cellDim[0] = (PRISM_FLOAT_PRECISION)val;
        std::cout << "Setting Z cell dimension to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setTileX_fromLineEdit(){
    int val = this->ui->lineEdit_tileX->text().toInt();
    if (val > 0){
        this->meta->tileX = (size_t)val;
        std::cout << "Setting tileX to " << val << " UCs" << std::endl;
        updateAlphaMax();
    }
    resetCalculation();
}

void PRISMMainWindow::setTileY_fromLineEdit(){
    int val = this->ui->lineEdit_tileY->text().toInt();
    if (val > 0){
        this->meta->tileY = (size_t)val;
        std::cout << "Setting tileY to " << val << " UCs" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setTileZ_fromLineEdit(){
    int val = this->ui->lineEdit_tileZ->text().toInt();
    if (val > 0){
        this->meta->tileZ = (size_t)val;
        std::cout << "Setting tileZ to " << val << " UCs" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setRandomSeed_fromLineEdit(){
    int val = this->ui->lineEdit_randomSeed->text().toInt();
    if (val >= 0){
        this->meta->random_seed = (size_t)val;
        std::cout << "Setting random seed to " << val << std::endl;
    }
    resetCalculation();
}


void PRISMMainWindow::setE0_fromLineEdit(){
    double val = this->ui->lineEdit_E0->text().toDouble();
    if (val > 0){
        meta->E0 = val * 1e3;
        std::cout << "Setting E0 to " << val << " keV" << std::endl;
        ui->lbl_lambda->setText(QString::fromUtf8("\u03BB = ") + QString::number(calculateLambda(*meta)) + QString::fromUtf8("\u212B"));
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_stepX_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeStepX->text().toDouble();
    if (val > 0){
        this->meta->probe_stepX = val;
        std::cout << "Setting probe_stepX to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_stepY_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeStepY->text().toDouble();
    if (val > 0){
        this->meta->probe_stepY = val;
        std::cout << "Setting probe_stepY to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_defocus_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeDefocus->text().toDouble();
    if (val >= 0){
        this->meta->probeDefocus = val * 1e-10;
        std::cout << "Setting probe defocus to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_C3_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_C3->text().toDouble();
    if (val >= 0){
        this->meta->C3 = val * 1e-10;
        std::cout << "Setting C3 to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_C5_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_C5->text().toDouble();
    if (val >= 0){
        this->meta->C5 = val * 1e-10;
        std::cout << "Setting C5 to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setDetector_angle_step_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_detectorAngle->text().toDouble();
    if (val > 0){
        this->meta->detector_angle_step = val / 1000;
        std::cout << "Setting detector angle step to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}


void PRISMMainWindow::setprobe_Xtilt_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeTiltX->text().toDouble();
    if (val > 0){
        this->meta->probeXtilt = val;
        std::cout << "Setting probe X tilt to " << val << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_Ytilt_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeTiltY->text().toDouble();
    if (val > 0){
        this->meta->probeYtilt = val;
        std::cout << "Setting probe Y tilt to " << val << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setBatchGPU_fromLineEdit(){
    int val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeTiltY->text().toInt();
    if (val > 0){
        this->meta->batch_size_GPU = val;
        std::cout << "Setting batch size (GPU) to " << val << std::endl;
    }
}

void PRISMMainWindow::setBatchCPU_fromLineEdit(){
    int val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_probeTiltY->text().toInt();
    if (val > 0){
        this->meta->batch_size_CPU = val;
        std::cout << "Setting batch size (CPU) to " << val << std::endl;
    }
}



void PRISMMainWindow::setscan_WindowXMin_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_scanWindowXMin->text().toDouble();
    val = std::min(this->meta->scanWindowXMax, std::max(val, (PRISM_FLOAT_PRECISION)0.0));
    if (val > 0){
        this->meta->scanWindowXMin = val;
        std::cout << "Setting scan window X min to " << val << std::endl;
    }
    ui->lineEdit_scanWindowXMin->setText(QString::number(val));
    resetCalculation();
}

void PRISMMainWindow::setscan_WindowXMax_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_scanWindowXMax->text().toDouble();
    val = std::max(this->meta->scanWindowXMin, std::min(val, (PRISM_FLOAT_PRECISION)1.0));
    if (val > 0){
        this->meta->scanWindowXMax = val;
        std::cout << "Setting scan window X max to " << val << std::endl;
    }
    ui->lineEdit_scanWindowXMax->setText(QString::number(val));
    resetCalculation();
}

void PRISMMainWindow::setscan_WindowYMin_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_scanWindowYMin->text().toDouble();
    val = std::min(this->meta->scanWindowYMax, std::max(val, (PRISM_FLOAT_PRECISION)0.0));
    if (val > 0){
        this->meta->scanWindowYMin = val;
        std::cout << "Setting scan window Y min to " << val << std::endl;
    }
    ui->lineEdit_scanWindowYMin->setText(QString::number(val));
    resetCalculation();
}

void PRISMMainWindow::setscan_WindowYMax_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_scanWindowYMax->text().toDouble();
    val = std::max(this->meta->scanWindowYMin, std::min(val, (PRISM_FLOAT_PRECISION)1.0));
    if (val > 0){
        this->meta->scanWindowYMax = val;
        std::cout << "Setting scan window Y max to " << val << std::endl;
    }
    ui->lineEdit_scanWindowYMax->setText(QString::number(val));
    resetCalculation();
}

void PRISMMainWindow::calculatePotential(){
    resetCalculation();
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();
    PotentialThread *worker = new PotentialThread(this, progressbar);
    worker->meta.toString();
    connect(worker, SIGNAL(finished()), this, SLOT(updatePotentialImage()));
    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
    worker->start();
}

void PRISMMainWindow::calculateSMatrix(){
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();
    SMatrixThread *worker = new SMatrixThread(this, progressbar);
    worker->meta.toString();
//    connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
    worker->start();
}

void PRISMMainWindow::calculateAll(){
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();

    if (meta->algorithm == PRISM::Algorithm::PRISM) {
	    FullPRISMCalcThread *worker = new FullPRISMCalcThread(this, progressbar);
	    connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
        connect(worker, SIGNAL(outputCalculated()), this, SLOT(updateOutputImage()));
        connect(worker, SIGNAL(outputCalculated()), this, SLOT(enableOutputWidgets()));
        connect(worker, SIGNAL(signalTitle(const QString)), progressbar, SLOT(setTitle(const QString)));
        connect(worker, SIGNAL(finished()), progressbar, SLOT(close()));
	    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
	    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
        std::cout <<"Starting Full PRISM Calculation" << std::endl;
        worker->meta.toString();
        worker->start();
    } else{
        FullMultisliceCalcThread *worker = new FullMultisliceCalcThread(this, progressbar);
        std::cout <<"Starting Full Multislice Calculation" << std::endl;
        worker->meta.toString();
        connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
        connect(worker, SIGNAL(outputCalculated()), this, SLOT(updateOutputImage()));
        connect(worker, SIGNAL(outputCalculated()), this, SLOT(enableOutputWidgets()));
        connect(worker, SIGNAL(finished()), progressbar, SLOT(close()));
        connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
        connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
        worker->start();
    }
}

void PRISMMainWindow::calculateProbe(){
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();
    PRISM_FLOAT_PRECISION X = (PRISM_FLOAT_PRECISION)ui->lineEdit_probeX->text().toDouble();
    PRISM_FLOAT_PRECISION Y = (PRISM_FLOAT_PRECISION)ui->lineEdit_probeY->text().toDouble();
    currently_calculated_X = X;
    currently_calculated_Y = Y;
    ProbeThread *worker = new ProbeThread(this, X, Y, progressbar, ui->checkBox_log->isChecked());
    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(worker, SIGNAL(finished()), this, SLOT(updatePotentialImage()));
//    connect(worker, SIGNAL(signal_pearsonReal(QString)), this, SLOT(update_pearsonReal(QString)));
//    connect(worker, SIGNAL(signal_pearsonK(QString)), this, SLOT(update_pearsonK(QString)));
    connect(worker, SIGNAL(signal_RReal(QString)), this, SLOT(update_RReal(QString)));
    connect(worker, SIGNAL(signal_RK(QString)), this, SLOT(update_RK(QString)));

//    connect(worker, SIGNAL(finished()), this, SLOT(updatePotentialDisplay()));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(close()));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
    connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
    connect(worker, SIGNAL(signalProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>)), this, SLOT(updateProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>)));
    connect(worker, SIGNAL(signalProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>)), this, SLOT(updateProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION>)));
    connect(worker, SIGNAL(signalProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>)), this, SLOT(updateProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>)));
    connect(worker, SIGNAL(signalProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>)), this, SLOT(updateProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION>)));
    connect(worker, SIGNAL(signalProbe_diffR(PRISM::Array2D<PRISM_FLOAT_PRECISION>, PRISM::Array2D<PRISM_FLOAT_PRECISION>)), this, SLOT(updateProbe_diffR(PRISM::Array2D<PRISM_FLOAT_PRECISION>, PRISM::Array2D<PRISM_FLOAT_PRECISION>)));
    connect(worker, SIGNAL(signalProbe_diffK(PRISM::Array2D<PRISM_FLOAT_PRECISION>, PRISM::Array2D<PRISM_FLOAT_PRECISION>)), this, SLOT(updateProbe_diffK(PRISM::Array2D<PRISM_FLOAT_PRECISION>, PRISM::Array2D<PRISM_FLOAT_PRECISION>)));
    worker->start();
}

void PRISMMainWindow::updatePotentialImage(){
    if (potentialIsReady()){
            {
            // create new empty image with appropriate dimensions
            potentialImage = QImage(pars.pot.get_dimj(), pars.pot.get_dimi(), QImage::Format_ARGB32);
            }

            // update sliders to match dimensions of potential, which also triggers a redraw of the image
            this->ui->slider_slicemin->setMinimum(1);
            this->ui->slider_slicemax->setMinimum(1);
            this->ui->slider_slicemin->setMaximum(pars.pot.get_dimk());
            this->ui->slider_slicemax->setMaximum(pars.pot.get_dimk());

            // I set the value to 0 and then to the correct value to ensure that the display update is triggered. A bit of a hack..
            this->ui->slider_slicemax->setValue(0);
            this->ui->slider_slicemax->setValue(pars.pot.get_dimk());
        }
}

void PRISMMainWindow::updatePotentialFloatImage(){
//    if (potentialReady){
    if (potentialImageExists){
//        QMutexLocker gatekeeper(&potentialLock);
        QMutexLocker gatekeeper(&dataLock);

        // integrate image into the float array, then convert to uchar
        size_t min_layer = this->ui->slider_slicemin->value();
        size_t max_layer = this->ui->slider_slicemax->value();
        potentialImage_float = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.pot.get_dimj(), pars.pot.get_dimi()}});
        for (auto k = min_layer; k <= max_layer; ++k){
            for (auto j = 0; j < pars.pot.get_dimj(); ++j){
                for (auto i = 0; i < pars.pot.get_dimi(); ++i){
                    potentialImage_float.at(j,i) += pars.pot.at(k - 1,j ,i);
                }
            }
        }

        // get max/min values for contrast setting
        auto minval = std::min_element(potentialImage_float.begin(),
                                                 potentialImage_float.end());
        auto maxval = std::max_element(potentialImage_float.begin(),
                                                 potentialImage_float.end());
        contrast_potentialMin = *minval;
        contrast_potentialMax = *maxval;
        ui->lineEdit_contrastPotMin->setText(QString::number(contrast_potentialMin));
        ui->lineEdit_contrastPotMax->setText(QString::number(contrast_potentialMax));
    }
    updatePotentialDisplay();
}

void PRISMMainWindow::updatePotentialDisplay(){
//    if (potentialReady){
    if (potentialImageExists){
//            QMutexLocker gatekeeper(&potentialLock);
            QMutexLocker gatekeeper(&dataLock);

            for (auto j = 0; j < pars.pot.get_dimj(); ++j){
                for (auto i = 0; i < pars.pot.get_dimi(); ++i){
                    uchar val = getUcharFromFloat(potentialImage_float.at(j,i),
                                                  contrast_potentialMin,
                                                  contrast_potentialMax);
                    potentialImage.setPixel(j, i, qRgba(val,val,val,255));
                }
            }

        QImage potentialImage_tmp = potentialImage.scaled(ui->lbl_image_potential->width(),
                               ui->lbl_image_potential->height(),
                               Qt::KeepAspectRatio);

        ui->lbl_image_potential->setPixmap(QPixmap::fromImage( potentialImage.scaled(ui->lbl_image_potential->width(),
                                                                                     ui->lbl_image_potential->height(),
                                                                                     Qt::KeepAspectRatio)));

        probeImage = potentialImage;
        PRISM_FLOAT_PRECISION xc, yc;
        std::cout <<"currently_calculated_X = " <<currently_calculated_X << std::endl;
        xc = currently_calculated_X / pars.pixelSize[1];
        yc = currently_calculated_Y / pars.pixelSize[0];

        long xc_im, yc_im;
        xc_im = (xc / pars.imageSize[1]) * probeImage.height();
        yc_im = (yc / pars.imageSize[0]) * probeImage.width();


        const long linehalfwidth = 1;
        const long linelength = 10;
//        std::cout << "test" << std::endl;
        for (auto ypix = -linehalfwidth; ypix <= linehalfwidth; ++ypix){
            for (auto x = -linelength; x <= linelength; ++x){
//                std::cout << "ypix + yc_im) % probeImage.width() = " << (ypix + yc_im) % probeImage.width()<< std::endl;
             probeImage.setPixel((probeImage.width()  + (ypix + yc_im) % probeImage.width()) % probeImage.width(),
                                 (probeImage.height() + (x + xc_im) % probeImage.height()) % probeImage.height(), qRgba(0, 255, 255, 100));
            }
        }

        for (auto xpix = -linehalfwidth; xpix <= linehalfwidth; ++xpix){
            for (auto y = -linelength; y <= linelength; ++y){
             probeImage.setPixel((probeImage.width()  + (y + yc_im) % probeImage.width()) % probeImage.width(),
                                 (probeImage.height() + (xpix + xc_im) % probeImage.height()) % probeImage.height(), qRgba(0, 255, 255, 100));
            }
        }

//        for (auto xpix = -linehalfwidth; xpix <= linehalfwidth; ++xpix){
//            for (auto y = 0; y < linelength; ++y){
//             probeImage.setPixel(xpix,y, qRgba(0, 255, 255, 100));
//            }
//        }

        ui->lbl_image_probeInteractive->setPixmap(QPixmap::fromImage( probeImage.scaled(ui->lbl_image_probeInteractive->width(),
                                                                                            ui->lbl_image_probeInteractive->height(),
                                                                                            Qt::KeepAspectRatio)));
    }
}

void PRISMMainWindow::updateProbeK_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION> arr){
    probeImage_pk = QImage(arr.get_dimj(), arr.get_dimi(), QImage::Format_ARGB32);
    auto contrast = std::minmax_element(arr.begin(), arr.end());
//    std::cout << "pK *contrast.first= " << *contrast.first<< std::endl;
//    std::cout << "pK *contrast.second= " << *contrast.first<< std::endl;
    for (auto j = 0; j < arr.get_dimj(); ++j){
        for (auto i = 0; i < arr.get_dimi(); ++i){
            uchar val = getUcharFromFloat(arr.at(j,i),*contrast.first, *contrast.second);
            probeImage_pk.setPixel(j, i, qRgba(val,val,val,255));
        }
    }
    ui->lbl_image_probe_pk->setPixmap(QPixmap::fromImage(probeImage_pk.scaled(ui->lbl_image_probe_pk->width(),
                                                                              ui->lbl_image_probe_pk->height(),
                                                                              Qt::KeepAspectRatio)));
    probeImage_pk_float = arr;
}
void PRISMMainWindow::updateProbeR_PRISM(PRISM::Array2D<PRISM_FLOAT_PRECISION> arr){
    probeImage_pr = QImage(arr.get_dimj(), arr.get_dimi(), QImage::Format_ARGB32);
    auto contrast = std::minmax_element(arr.begin(), arr.end());
//    std::cout << "pReal *contrast.first= " << *contrast.first<< std::endl;
//    std::cout << "pReal *contrast.second= " << *contrast.first<< std::endl;
    for (auto j = 0; j < arr.get_dimj(); ++j){
        for (auto i = 0; i < arr.get_dimi(); ++i){
            uchar val = getUcharFromFloat(arr.at(j,i),*contrast.first, *contrast.second);
            probeImage_pr.setPixel(j, i, qRgba(val,val,val,255));
        }
    }
    ui->lbl_image_probe_pr->setPixmap(QPixmap::fromImage(probeImage_pr.scaled(ui->lbl_image_probe_pr->width(),
                                                                              ui->lbl_image_probe_pr->height(),
                                                                              Qt::KeepAspectRatio)));
    probeImage_pr_float = arr;
}
void PRISMMainWindow::updateProbeK_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION> arr){
    probeImage_mk = QImage(arr.get_dimj(), arr.get_dimi(), QImage::Format_ARGB32);
    auto contrast = std::minmax_element(arr.begin(), arr.end());
//    std::cout << "mK *contrast.first= " << *contrast.first<< std::endl;
//    std::cout << "mK *contrast.second= " << *contrast.first<< std::endl;
    for (auto j = 0; j < arr.get_dimj(); ++j){
        for (auto i = 0; i < arr.get_dimi(); ++i){
            uchar val = getUcharFromFloat(arr.at(j,i),*contrast.first, *contrast.second);
            probeImage_mk.setPixel(j, i, qRgba(val,val,val,255));
        }
    }
    ui->lbl_image_probe_mk->setPixmap(QPixmap::fromImage(probeImage_mk.scaled(ui->lbl_image_probe_mk->width(),
                                                                              ui->lbl_image_probe_mk->height(),
                                                                              Qt::KeepAspectRatio)));
    probeImage_mk_float = arr;
}
void PRISMMainWindow::updateProbeR_Multislice(PRISM::Array2D<PRISM_FLOAT_PRECISION> arr){
    probeImage_mr = QImage(arr.get_dimj(), arr.get_dimi(), QImage::Format_ARGB32);
    auto contrast = std::minmax_element(arr.begin(), arr.end());
//    std::cout << "mReal *contrast.first= " << *contrast.first<< std::endl;
//    std::cout << "mReal *contrast.second= " << *contrast.first<< std::endl;
    for (auto j = 0; j < arr.get_dimj(); ++j){
        for (auto i = 0; i < arr.get_dimi(); ++i){
            uchar val = getUcharFromFloat(arr.at(j,i),*contrast.first, *contrast.second);
            probeImage_mr.setPixel(j, i, qRgba(val,val,val,255));
        }
    }
    ui->lbl_image_probe_mr->setPixmap(QPixmap::fromImage(probeImage_mr.scaled(ui->lbl_image_probe_mr->width(),
                                                                              ui->lbl_image_probe_mr->height(),
                                                                              Qt::KeepAspectRatio)));
    probeImage_mr_float = arr;
}

void PRISMMainWindow::updateProbe_diffR(PRISM::Array2D<PRISM_FLOAT_PRECISION> arr, PRISM::Array2D<PRISM_FLOAT_PRECISION> arr_contrast){
    probeImage_diffr = QImage(arr.get_dimj(), arr.get_dimi(), QImage::Format_ARGB32);
    auto contrast = std::minmax_element(arr_contrast.begin(), arr_contrast.end());
//    std::cout << "diffreal *contrast.first= " << *contrast.first<< std::endl;
//    std::cout << "diffreal *contrast.second= " << *contrast.first<< std::endl;
    for (auto j = 0; j < arr.get_dimj(); ++j){
        for (auto i = 0; i < arr.get_dimi(); ++i){
            uchar val = getUcharFromFloat(arr.at(j,i),*contrast.first, *contrast.second);
            probeImage_diffr.setPixel(j, i, qRgba(val,val,val,255));
        }
    }
    ui->lbl_image_probeDifferenceR->setPixmap(QPixmap::fromImage(probeImage_diffr.scaled(ui->lbl_image_probeDifferenceR->width(),
                                                                                         ui->lbl_image_probeDifferenceR->height(),
                                                                                         Qt::KeepAspectRatio)));
    probeImage_diffr_float = arr;
}
void PRISMMainWindow::updateProbe_diffK(PRISM::Array2D<PRISM_FLOAT_PRECISION> arr, PRISM::Array2D<PRISM_FLOAT_PRECISION> arr_contrast){
    probeImage_diffk = QImage(arr.get_dimj(), arr.get_dimi(), QImage::Format_ARGB32);
    auto contrast = std::minmax_element(arr_contrast.begin(), arr_contrast.end());
//    std::cout << "diffk *contrast.first= " << *contrast.first<< std::endl;
//    std::cout << "diffk *contrast.second= " << *contrast.first<< std::endl;

    for (auto j = 0; j < arr.get_dimj(); ++j){
        for (auto i = 0; i < arr.get_dimi(); ++i){
            uchar val = getUcharFromFloat(arr.at(j,i),*contrast.first, *contrast.second);
            probeImage_diffk.setPixel(j, i, qRgba(val,val,val,255));
        }
    }
    ui->lbl_image_probeDifferenceK->setPixmap(QPixmap::fromImage(probeImage_diffk.scaled(ui->lbl_image_probeDifferenceK->width(),
                                                                                         ui->lbl_image_probeDifferenceK->height(),
                                                                                         Qt::KeepAspectRatio)));
    probeImage_diffk_float = arr;
}

void PRISMMainWindow::updateOutputImage(){
//    if (outputReady){
    if (outputReady){
	    std::cout << "updateOutputImage " << std::endl;
            {
//            QMutexLocker gatekeeper(&outputLock);
            QMutexLocker gatekeeper(&dataLock);

            // create new empty image with appropriate dimensions
            outputImage = QImage(pars.output.get_dimk(), pars.output.get_dimj(), QImage::Format_ARGB32);
            }
            // update sliders to match dimensions of output, which also triggers a redraw of the image
            this->ui->slider_angmin->setMinimum(0);
            this->ui->slider_angmax->setMinimum(0);
            this->ui->slider_angmin->setMaximum(detectorAngles.size() - 1);
            this->ui->slider_angmax->setMaximum(detectorAngles.size() - 1);
            this->ui->slider_angmax->setValue(detectorAngles.size() - 1);
	        this->ui->lineEdit_angmin->setText(QString::number(detectorAngles[0]));
	        this->ui->lineEdit_angmax->setText(QString::number(detectorAngles[detectorAngles.size() - 1]));
        }
    updateOutputFloatImage();
}

void PRISMMainWindow::updateOutputFloatImage(){
//    if (outputReady){
      if (outputImageExists){
        std::cout << "updateOutputFloatImage " << std::endl;

//        QMutexLocker gatekeeper(&outputLock);
        QMutexLocker gatekeeper(&dataLock);

        // integrate image into the float array, then convert to uchar
        size_t min_layer = this->ui->slider_angmin->value();
        size_t max_layer = this->ui->slider_angmax->value();
        std::cout << "min_layer = " << min_layer << std::endl;
        std::cout << "max_layer = " << max_layer << std::endl;
        outputImage_float = PRISM::zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.output.get_dimk(), pars.output.get_dimj()}});
//        std::cout << "outputImage_float.get_dimj() = " << outputImage_float.get_dimj() << std::endl;
//         std::cout << "outputImage_float.get_dimi() = " << outputImage_float.get_dimi() << std::endl;
        for (auto j = 0; j < pars.output.get_dimk(); ++j){
            for (auto i = 0; i < pars.output.get_dimj(); ++i){
                 for (auto k = min_layer; k <= max_layer; ++k){
                    outputImage_float.at(j,i) += pars.output.at(j, i, k);
                }
            }
        }

        // get max/min values for contrast setting
        auto minval = std::min_element(outputImage_float.begin(),
                                       outputImage_float.end());
        auto maxval = std::max_element(outputImage_float.begin(),
                                       outputImage_float.end());
        contrast_outputMin = *minval;
        contrast_outputMax = *maxval;
        ui->contrast_outputMin->setText(QString::number(contrast_outputMin));
        ui->contrast_outputMax->setText(QString::number(contrast_outputMax));
    }
    updateOutputDisplay();
}

void PRISMMainWindow::updateOutputDisplay(){
//    if (outputReady){
    if (outputImageExists){
//            QMutexLocker gatekeeper(&outputLock);
        QMutexLocker gatekeeper(&dataLock);
            for (auto j = 0; j < pars.output.get_dimk(); ++j){
                for (auto i = 0; i < pars.output.get_dimj(); ++i){
                    uchar val = getUcharFromFloat(outputImage_float.at(j,i),
                                                  contrast_outputMin,
                                                  contrast_outputMax);
                    outputImage.setPixel(j, i, qRgba(val,val,val,255));
                }
            }

        QImage outputImage_tmp = outputImage.scaled(ui->lbl_image_output->width(),
                                                    ui->lbl_image_output->height(),
                                                    Qt::KeepAspectRatio);

        ui->lbl_image_output->setPixmap(QPixmap::fromImage( outputImage.scaled(ui->lbl_image_output->width(),
                                                                               ui->lbl_image_output->height(),
                                                                               Qt::KeepAspectRatio)));
    }
}

void PRISMMainWindow::updateSliders_fromLineEdits(){
    this->ui->slider_slicemin->setValue(std::min(this->ui->lineEdit_slicemin->text().toInt(),
                                                 this->ui->slider_slicemax->value()));
    this->ui->slider_slicemax->setValue(std::max(this->ui->lineEdit_slicemax->text().toInt(),
                                                 this->ui->slider_slicemin->value()));
}
void PRISMMainWindow::updateSliders_fromLineEdits_ang(){
//    if (outputReady){
    if (outputImageExists){
        PRISM_FLOAT_PRECISION minval = ( (PRISM_FLOAT_PRECISION)this->ui->lineEdit_angmin->text().toDouble() - detectorAngles[0]) /
        (detectorAngles[1]-detectorAngles[0]);
        PRISM_FLOAT_PRECISION maxval = ( (PRISM_FLOAT_PRECISION)this->ui->lineEdit_angmax->text().toDouble() - detectorAngles[0]) /
        (detectorAngles[1]-detectorAngles[0]);
        this->ui->slider_angmin->setValue(std::min( (int)std::round(minval),
                                                     this->ui->slider_angmax->value()));
        this->ui->slider_angmax->setValue(std::max( (int)std::round(maxval),
                                                     this->ui->slider_angmin->value()));
        updateSlider_lineEdits_max_ang(this->ui->slider_angmax->value());
        updateSlider_lineEdits_min_ang(this->ui->slider_angmin->value());
    }
}

void PRISMMainWindow::updateSlider_lineEdits_min(int val){
    if (val <= this->ui->slider_slicemax->value()){
        this->ui->lineEdit_slicemin->setText(QString::number(val));
    } else {
        this->ui->slider_slicemin->setValue(this->ui->slider_slicemax->value());
    }
}

void PRISMMainWindow::updateSlider_lineEdits_max(int val){
    if (val >= this->ui->slider_slicemin->value()){
        this->ui->lineEdit_slicemax->setText(QString::number(val));
    } else {
        this->ui->slider_slicemax->setValue(this->ui->slider_slicemin->value());
    }
}

void PRISMMainWindow::updateSlider_lineEdits_max_ang(int val){
//    if (outputReady){
    QMutexLocker gatekeeper(&dataLock);
    if (outputImageExists){
        if (val >= this->ui->slider_angmin->value()){
            double scaled_val = detectorAngles[0] + val * (detectorAngles[1] - detectorAngles[0]);
            std::cout << "val = " << val << std::endl;
            std::cout << "scaled_val = " << scaled_val << std::endl;
            this->ui->lineEdit_angmax->setText(QString::number(scaled_val));
        } else {
            this->ui->slider_angmax->setValue(this->ui->slider_angmin->value());
        }
    }
}

void PRISMMainWindow::updateSlider_lineEdits_min_ang(int val){
//    if (outputReady){
    QMutexLocker gatekeeper(&dataLock);
    if (outputImageExists){
        if (val <= this->ui->slider_angmax->value()){
            double scaled_val = detectorAngles[0] + val * (detectorAngles[1] - detectorAngles[0]);
            this->ui->lineEdit_angmin->setText(QString::number(scaled_val));
        } else {
            this->ui->slider_angmin->setValue(this->ui->slider_angmax->value());
        }
    }
}


void PRISMMainWindow::updateContrastPotMin(){
    contrast_potentialMin = (PRISM_FLOAT_PRECISION)ui->lineEdit_contrastPotMin->text().toDouble();
    updatePotentialDisplay();
}
void PRISMMainWindow::updateContrastPotMax(){
    contrast_potentialMax = (PRISM_FLOAT_PRECISION)ui->lineEdit_contrastPotMax->text().toDouble();
    updatePotentialDisplay();

}

void PRISMMainWindow::updateContrastAngMin(){
    contrast_outputMin = (PRISM_FLOAT_PRECISION)ui->contrast_outputMin->text().toDouble();
    updateOutputDisplay();
}
void PRISMMainWindow::updateContrastAngMax(){
    contrast_outputMax = (PRISM_FLOAT_PRECISION)ui->contrast_outputMax->text().toDouble();
    updateOutputDisplay();
}

void PRISMMainWindow::updateAlphaMax(){
    using namespace PRISM;
    PRISM_FLOAT_PRECISION f_x = 4 * meta->interpolationFactorX;
    PRISM_FLOAT_PRECISION f_y = 4 * meta->interpolationFactorY;
    Array1D<size_t> imageSize({{(size_t)meta->cellDim[1] * meta->tileY, (size_t)meta->cellDim[2] * meta->tileX}}, {{2}});
    imageSize[0] = (size_t)std::max((PRISM_FLOAT_PRECISION)4.0,  (PRISM_FLOAT_PRECISION)(f_y * round(((PRISM_FLOAT_PRECISION)imageSize[0]) / meta->realspace_pixelSize / f_y)));
    imageSize[1] = (size_t)std::max((PRISM_FLOAT_PRECISION)4.0,  (PRISM_FLOAT_PRECISION)(f_x * round(((PRISM_FLOAT_PRECISION)imageSize[1]) / meta->realspace_pixelSize / f_x)));

    Array1D<PRISM_FLOAT_PRECISION> qx = makeFourierCoords(imageSize[1], meta->realspace_pixelSize);
    Array1D<PRISM_FLOAT_PRECISION> qy = makeFourierCoords(imageSize[0], meta->realspace_pixelSize);
    std::pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh = meshgrid(qy, qx);

    PRISM_FLOAT_PRECISION qMax;
    Array2D<PRISM_FLOAT_PRECISION> qxa;
    Array2D<PRISM_FLOAT_PRECISION> qya;

    qya = mesh.first;
    qxa = mesh.second;
    // get qMax
    qMax = 0;
    {
        PRISM_FLOAT_PRECISION qx_max = 0;
        PRISM_FLOAT_PRECISION qy_max = 0;
        for (auto i = 0; i < qx.size(); ++i) {
            qx_max = (std::abs(qx[i]) > qx_max) ? std::abs(qx[i]) : qx_max;
            qy_max = (std::abs(qy[i]) > qy_max) ? std::abs(qy[i]) : qy_max;
        }

        long long ncx = (size_t) floor((PRISM_FLOAT_PRECISION) imageSize[1] / 2);
        PRISM_FLOAT_PRECISION dpx = 1 / (imageSize[1] * meta->realspace_pixelSize);
        long long ncy = (size_t) floor((PRISM_FLOAT_PRECISION) imageSize[0] / 2);
        PRISM_FLOAT_PRECISION dpy = 1 / (imageSize[0] * meta->realspace_pixelSize);
        qMax = std::min(qx_max, qy_max) / 2;
    }


//    long long ncx = (size_t) floor((PRISM_FLOAT_PRECISION) imageSize[1] / 2);
//    PRISM_FLOAT_PRECISION dpx = 1 / (imageSize[1] * meta->realspace_pixelSize);
//    long long ncy = (size_t) floor((PRISM_FLOAT_PRECISION) imageSize[0] / 2);
//    PRISM_FLOAT_PRECISION dpy = 1 / (imageSize[0] * meta->realspace_pixelSize);
//    PRISM_FLOAT_PRECISION qMax = std::min(dpx*(ncx - 1), dpy*(ncy - 1)) / 2;


    PRISM_FLOAT_PRECISION alphaMax = qMax * calculateLambda(*meta);
//    ui->lbl_alphaMax->setText("test");
    ui->lbl_alphaMax->setText(QString::fromUtf8("\u03B1 max = ") + QString::number(alphaMax));
}

void PRISMMainWindow::saveCurrentOutputImage(){
    if (outputImageExists){
        std::cout << "saving\n";
//            QMutexLocker gatekeeper(&outputLock);
        QMutexLocker gatekeeper(&dataLock);
        outputImage_float.toMRC_f(ui->lineEdit_saveOutputImage->text().toStdString().c_str());
    }
}

void PRISMMainWindow::toggleStreamingMode(){
    meta->transfer_mode = ui->checkBox_streamdata->isChecked() ? PRISM::StreamingMode::Stream : PRISM::StreamingMode::SingleXfer;
}

void PRISMMainWindow::newRandomSeed(){
    PRISM_FLOAT_PRECISION val = rand() % 100000;
    ui->lineEdit_randomSeed->setText(QString::number(val));
    meta->random_seed = val;
}

void PRISMMainWindow::toggle3DOutput(){
    meta->save3DOutput = ui->checkBox_3D->isChecked();
}

void PRISMMainWindow::toggle4DOutput(){
    meta->save4DOutput = ui->checkBox_4D->isChecked();
}

void PRISMMainWindow::toggleThermalEffects(){
    meta->include_thermal_effects = ui->checkBox_thermalEffects->isChecked();
    resetCalculation();
}

void PRISMMainWindow::update_pearsonReal(QString str){
    ui->lbl_pears_real->setText(str);
}

void PRISMMainWindow::update_pearsonK(QString str){
    ui->lbl_pears_k->setText(str);
}

void PRISMMainWindow::update_RReal(QString str){
    ui->lbl_R_real->setText(str);
}

void PRISMMainWindow::update_RK(QString str){
    ui->lbl_R_k->setText(str);
}


void PRISMMainWindow::toggleSaveProjectedPotential(){
    this->saveProjectedPotential = ui->checkBox_saveProjectedPotential->isChecked() ? true:false;
}

bool PRISMMainWindow::potentialIsReady(){
    QMutexLocker gatekeeper(&dataLock);
    return potentialReady;
}
bool PRISMMainWindow::SMatrixIsReady(){
    QMutexLocker gatekeeper(&dataLock);
    return ScompactReady;
}
bool PRISMMainWindow::OutputIsReady(){
    QMutexLocker gatekeeper(&dataLock);
    return outputReady;
}

void PRISMMainWindow::enableOutputWidgets(){
    ui->slider_angmax->setEnabled(true);
    ui->slider_angmin->setEnabled(true);
    ui->contrast_outputMax->setEnabled(true);
    ui->contrast_outputMin->setEnabled(true);
    ui->lineEdit_angmin->setEnabled(true);
    ui->lineEdit_angmax->setEnabled(true);
}
void PRISMMainWindow::resetCalculation(){
//	QMutexLocker gatekeeper(&calculationLock);
    QMutexLocker gatekeeper(&dataLock);
    std::cout << "Resetting Calculation" << std::endl;
    potentialReady  = false;
    ScompactReady   = false;
    outputReady     = false;
    probeSetupReady = false;
}

void PRISMMainWindow::redrawImages(){
    updatePotentialDisplay();
//    ui->lbl_image_potential->setPixmap(QPixmap::fromImage(potentialImage.scaled(ui->lbl_image_potential->width(),
//                                                                                ui->lbl_image_potential->height(),
//                                                                                Qt::KeepAspectRatio)));

    ui->lbl_image_probeInteractive->setPixmap(QPixmap::fromImage(probeImage.scaled(ui->lbl_image_probeInteractive->width(),
                                                                                   ui->lbl_image_probeInteractive->height(),
                                                                                   Qt::KeepAspectRatio)));

    ui->lbl_image_output->setPixmap(QPixmap::fromImage(outputImage.scaled(ui->lbl_image_output->width(),
                                                                          ui->lbl_image_output->height(),
                                                                          Qt::KeepAspectRatio)));

    ui->lbl_image_probe_pk->setPixmap(QPixmap::fromImage(probeImage_pk.scaled(ui->lbl_image_probe_pk->width(),
                                                                              ui->lbl_image_probe_pk->height(),
                                                                              Qt::KeepAspectRatio)));

    ui->lbl_image_probe_pr->setPixmap(QPixmap::fromImage(probeImage_pr.scaled(ui->lbl_image_probe_pr->width(),
                                                                              ui->lbl_image_probe_pr->height(),
                                                                              Qt::KeepAspectRatio)));

    ui->lbl_image_probe_mk->setPixmap(QPixmap::fromImage(probeImage_mk.scaled(ui->lbl_image_probe_mk->width(),
                                                                              ui->lbl_image_probe_mk->height(),
                                                                              Qt::KeepAspectRatio)));

    ui->lbl_image_probe_mr->setPixmap(QPixmap::fromImage(probeImage_mr.scaled(ui->lbl_image_probe_mr->width(),
                                                                              ui->lbl_image_probe_mr->height(),
                                                                              Qt::KeepAspectRatio)));
    ui->lbl_image_probeDifferenceK->setPixmap(QPixmap::fromImage(probeImage_diffk.scaled(ui->lbl_image_probeDifferenceK->width(),
                                                                                         ui->lbl_image_probeDifferenceK->height(),
                                                                                         Qt::KeepAspectRatio)));

    ui->lbl_image_probeDifferenceR->setPixmap(QPixmap::fromImage(probeImage_diffr.scaled(ui->lbl_image_probeDifferenceR->width(),
                                                                                         ui->lbl_image_probeDifferenceR->height(),
                                                                                         Qt::KeepAspectRatio)));
    ui->lbl_image_output->setPixmap(QPixmap::fromImage(outputImage.scaled(ui->lbl_image_output->width(),
                                                                          ui->lbl_image_output->height(),
                                                                          Qt::KeepAspectRatio)));

//    std::cout << "ui->lbl_image_output->width() = " << ui->lbl_image_output->width() << std::endl;
//    std::cout << "ui->lbl_image_output->height() = " << ui->lbl_image_output->height() << std::endl;

}


void PRISMMainWindow::resizeEvent(QResizeEvent* event)
{
   QMainWindow::resizeEvent(event);
   redrawImages();
}


PRISMMainWindow::~PRISMMainWindow()
{
    delete ui;
	delete meta;
}

unsigned char getUcharFromFloat(PRISM_FLOAT_PRECISION val,
                                PRISM_FLOAT_PRECISION contrast_low,
                                PRISM_FLOAT_PRECISION contrast_high){

    if (val < contrast_low)  return 0;
    if (val > contrast_high) return 255;
//    return (unsigned char)( val / contrast_high * 255);
    return (unsigned char)( (val - contrast_low) / (contrast_high - contrast_low) * 255);

}

PRISM_FLOAT_PRECISION calculateLambda(PRISM::Metadata<PRISM_FLOAT_PRECISION> meta){
	constexpr double m = 9.109383e-31;
	constexpr double e = 1.602177e-19;
	constexpr double c = 299792458;
	constexpr double h = 6.62607e-34;
	return (PRISM_FLOAT_PRECISION)(h / sqrt(2 * m * e * meta.E0) / sqrt(1 + e * meta.E0 / 2 / m / c / c) * 1e10);
}
