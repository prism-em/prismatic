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

#include <QFileDialog>
#include <QMessageBox>
#include <QPainter>
#include "prismmainwindow.h"
#include "ui_prismmainwindow.h"
#include <fstream>
#include <iostream>
#include <utility>
#include <array>
#include "configure.h"
#include "prism_qthreads.h"
#include "prism_progressbar.h"
#include "saveatomiccoordinatesdialog.h"
#include "utility.h"
#include "atom.h"
#include "parseInput.h"
#include "params.h"
#include "aberration.h"
#include "probe.h"
#include <cstdio>
#include "QMessageBox"
#include <stdio.h>
#include <QApplication>
#include <QFile>
#include <QTextStream>
#include <QPropertyAnimation>
//#include <unistd.h>

bool validateFilename(const std::string str){
    std::ifstream f(str);
    return f.good();
}

bool validateWriteFilename(const std::string str){
    std::cout << "Validating file " << str << " for writing" << std::endl;
    std::ofstream f(str);
    if (f.good())std::cout <<"file " << str << " good" << std::endl;
    return f.good();
}

std::string get_default_parameter_filename() {

#ifdef _WIN32
	char* appdata = getenv("APPDATA");
	return std::string(appdata) + "\\prismatic_gui_params.txt";
#else
	char* appdata = getenv("HOME");
	return std::string(appdata) + "/prismatic_gui_params.txt";
#endif //_WIN32
}

PRISMATIC_FLOAT_PRECISION calculateLambda(Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta);

PRISMMainWindow::PRISMMainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::PRISMMainWindow),
    potentialReady(false),
    ScompactReady(false),
    outputReady(false),
    //saveProjectedPotential(false),
    probeSetupReady(false),
    potentialArrayExists(false),
    outputArrayExists(false),
    interpYSet(false),
    pixelSizeYSet(false),
    probeStepYSet(false),
    probeTiltYSet(false),
    minWindowYSet(false),
    maxWindowYSet(false),
    potentialImage(QImage()),
    currently_calculated_X(0.0),
    currently_calculated_Y(0.0),
    pixelSize({ 1,1 }),
    colormapper(Prismatic::Colormapper(Prismatic::VioletFireColormap))
{
    qRegisterMetaType<Prismatic::Array2D< PRISMATIC_FLOAT_PRECISION> >("Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>");
    qRegisterMetaType<Prismatic::Array3D< PRISMATIC_FLOAT_PRECISION> >("Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION>");

    // build Qt generated interface
    ui->setupUi(this);

    // set window title
    setWindowTitle("Prismatic (no atomic coordinate file selected)");

   /* ui->box_samplesettings->setStyleSheet("QGroupBox { \
                                          border: 1px solid gray;\
                                          border-radius: 9px;\
                                          margin-top: 0.5em;\
                                      }  QGroupBox::title {\
                                         subcontrol-origin: margin;\
                                         left: 150px;\
                                         padding: 0 3px 0 3px;\
                                         }");

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


    ui->box_simulationsettings->setStyleSheet("QGroupBox { \
                                        border: 1px solid gray;\
                                        border-radius: 9px;\
                                        margin-top: 0.5em;\
                                    }  QGroupBox::title {\
                                       subcontrol-origin: margin;\
                                       left: 150px;\
                                       padding: 0 3px 0px 3px;\
                                       }");*/


    //potentialImage.load(":/images/prism.png");
    potentialImage.load(":/images/prismatic-logo.png");
    probeImage.load(":/images/prismatic-logo.png");
    outputImage.load(":/images/prismatic-logo.png");
    outputImage_HRTEM.load(":/images/prismatic-logo.png");
    probeImage_pr.load(":/images/airy.png");
    probeImage_pk.load(":/images/airy.png");
    probeImage_mr.load(":/images/airy.png");
    probeImage_mk.load(":/images/airy.png");
    probeImage_diffk.load(":/images/airy.png");
    probeImage_diffr.load(":/images/airy.png");
    ui->lbl_image_potential->setPixmap(QPixmap::fromImage(potentialImage.scaled(540,
        420,
        Qt::KeepAspectRatio)));
    ui->lbl_image_potential_2->setPixmap(QPixmap::fromImage(potentialImage.scaled(540,
        420,
        Qt::KeepAspectRatio)));

    //    probeImage.load(":/images/probe.png");
    //    outputImage.load(":/images/output.png");
    redrawImages();


    // set initially displayed values based on the default parameters
    this->meta = new Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION>;
    updateDisplay();



    if (validateFilename(get_default_parameter_filename())) {
        readParams(get_default_parameter_filename());
    }

    // connect signals and slots
    connect(this->ui->btn_loadParams, SIGNAL(pressed()), this, SLOT(selectParameterFile()));
    connect(this->ui->actionLoad_Parameters, SIGNAL(triggered()), this, SLOT(selectParameterFile()));
    connect(this->ui->btn_saveParams, SIGNAL(pressed()), this, SLOT(writeParameterFile()));
    connect(this->ui->actionSave_Parameters, SIGNAL(triggered()), this, SLOT(writeParameterFile()));
    connect(this->ui->lineEdit_interpFactor_x, SIGNAL(textEdited(QString)), this, SLOT(setInterpolationFactorX()));
    connect(this->ui->lineEdit_interpFactor_y, SIGNAL(textEdited(QString)), this, SLOT(setInterpolationFactorY()));
    connect(this->ui->lineEdit_outputfile, SIGNAL(textEdited(QString)), this, SLOT(setFilenameOutput_fromLineEdit()));
    connect(this->ui->btn_atomsfile_browse, SIGNAL(pressed()), this, SLOT(setFilenameAtoms_fromDialog()));
    connect(this->ui->actionLoad_Coordinates, SIGNAL(triggered()), this, SLOT(setFilenameAtoms_fromDialog()));
    connect(this->ui->btn_saveCoordinates, SIGNAL(pressed()), this, SLOT(openSaveAtomsDialog()));
    connect(this->ui->actionSave_Coordinates, SIGNAL(triggered()), this, SLOT(openSaveAtomsDialog()));
    connect(this->ui->btn_loadAbFile, SIGNAL(pressed()), this, SLOT(readAberrationFile()));
    connect(this->ui->btn_loadProbePos, SIGNAL(pressed()), this, SLOT(readProbeFile()));
    connect(this->ui->spinBox_numGPUs, SIGNAL(valueChanged(int)), this, SLOT(setNumGPUs(const int&)));
    connect(this->ui->spinBox_numThreads, SIGNAL(valueChanged(int)), this, SLOT(setNumThreads(const int&)));
    connect(this->ui->spinBox_numFP, SIGNAL(valueChanged(int)), this, SLOT(setNumFP(const int&)));
    connect(this->ui->spinBox_numNS, SIGNAL(valueChanged(int)), this, SLOT(setNumNS(const int&)));
    connect(this->ui->spinBox_zSampling, SIGNAL(valueChanged(int)), this, SLOT(setzSampling(const int&)));
    connect(this->ui->spinBox_numStreams, SIGNAL(valueChanged(int)), this, SLOT(setNumStreams(const int&)));
    connect(this->ui->lineEdit_probeSemiangle, SIGNAL(textEdited(QString)), this, SLOT(setprobeSemiangle_fromLineEdit()));
    connect(this->ui->lineEdit_zStart, SIGNAL(textEdited(QString)), this, SLOT(setzStart_fromLineEdit()));
    // connect(this->ui->lineEdit_alphaBeamMax, SIGNAL(textEdited(QString)), this, SLOT(setalphaBeamMax_fromLineEdit()));
    connect(this->ui->lineEdit_2D_inner, SIGNAL(textEdited(QString)), this, SLOT(set2D_innerAngle_fromLineEdit()));
    connect(this->ui->lineEdit_2D_outer, SIGNAL(textEdited(QString)), this, SLOT(set2D_outerAngle_fromLineEdit()));
    connect(this->ui->lineEdit_pixelSizeX, SIGNAL(textEdited(QString)), this, SLOT(setPixelSizeX_fromLineEdit()));
    connect(this->ui->lineEdit_pixelSizeY, SIGNAL(textEdited(QString)), this, SLOT(setPixelSizeY_fromLineEdit()));
    connect(this->ui->lineEdit_batchCPU, SIGNAL(textEdited(QString)), this, SLOT(setBatchCPU_fromLineEdit()));
    connect(this->ui->lineEdit_batchGPU, SIGNAL(textEdited(QString)), this, SLOT(setBatchGPU_fromLineEdit()));
    connect(this->ui->lineEdit_potbound, SIGNAL(textEdited(QString)), this, SLOT(setPotBound_fromLineEdit()));
    connect(this->ui->lineEdit_sliceThickness, SIGNAL(textEdited(QString)), this, SLOT(setSliceThickness_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimX, SIGNAL(textEdited(QString)), this, SLOT(setCellDimX_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimY, SIGNAL(textEdited(QString)), this, SLOT(setCellDimY_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimZ, SIGNAL(textEdited(QString)), this, SLOT(setCellDimZ_fromLineEdit()));
    connect(this->ui->lineEdit_tileX, SIGNAL(textEdited(QString)), this, SLOT(setTileX_fromLineEdit()));
    connect(this->ui->lineEdit_tileY, SIGNAL(textEdited(QString)), this, SLOT(setTileY_fromLineEdit()));
    connect(this->ui->lineEdit_tileZ, SIGNAL(textEdited(QString)), this, SLOT(setTileZ_fromLineEdit()));
    connect(this->ui->lineEdit_randomSeed, SIGNAL(textEdited(QString)), this, SLOT(setRandomSeed_fromLineEdit()));
    connect(this->ui->lineEdit_probeDefocus, SIGNAL(textEdited(QString)), this, SLOT(setprobe_defocus_fromLineEdit()));
    connect(this->ui->lineEdit_dfr_min, SIGNAL(textEdited(QString)), this, SLOT(set_dfr_min_fromLineEdit()));
    connect(this->ui->lineEdit_dfr_max, SIGNAL(textEdited(QString)), this, SLOT(set_dfr_max_fromLineEdit()));
    connect(this->ui->lineEdit_dfr_step, SIGNAL(textEdited(QString)), this, SLOT(set_dfr_step_fromLineEdit()));
    connect(this->ui->lineEdit_dfs, SIGNAL(textEdited(QString)), this, SLOT(set_dfs_fromLineEdit()));
    connect(this->ui->lineEdit_C3, SIGNAL(textEdited(QString)), this, SLOT(setprobe_C3_fromLineEdit()));
    connect(this->ui->lineEdit_C5, SIGNAL(textEdited(QString)), this, SLOT(setprobe_C5_fromLineEdit()));
    connect(this->ui->lineEdit_detectorAngle, SIGNAL(textEdited(QString)), this, SLOT(setdetectorAngleStep_fromLineEdit()));
    connect(this->ui->lineEdit_probeTiltX, SIGNAL(textEdited(QString)), this, SLOT(setprobe_Xtilt_fromLineEdit()));
    connect(this->ui->lineEdit_probeTiltY, SIGNAL(textEdited(QString)), this, SLOT(setprobe_Ytilt_fromLineEdit()));
    connect(this->ui->lineEdit_xtt_min, SIGNAL(textEdited(QString)), this, SLOT(setxtt_min_fromLineEdit())); 
    connect(this->ui->lineEdit_xtt_max, SIGNAL(textEdited(QString)), this, SLOT(setxtt_max_fromLineEdit())); 
    connect(this->ui->lineEdit_xtt_step, SIGNAL(textEdited(QString)), this, SLOT(setxtt_step_fromLineEdit())); 
    connect(this->ui->lineEdit_ytt_min, SIGNAL(textEdited(QString)), this, SLOT(setytt_min_fromLineEdit())); 
    connect(this->ui->lineEdit_ytt_max, SIGNAL(textEdited(QString)), this, SLOT(setytt_max_fromLineEdit())); 
    connect(this->ui->lineEdit_ytt_step, SIGNAL(textEdited(QString)), this, SLOT(setytt_step_fromLineEdit())); 
    connect(this->ui->lineEdit_rtt_min, SIGNAL(textEdited(QString)), this, SLOT(setrtt_min_fromLineEdit())); 
    connect(this->ui->lineEdit_rtt_max, SIGNAL(textEdited(QString)), this, SLOT(setrtt_max_fromLineEdit())); 
    connect(this->ui->lineEdit_xtilt_offset, SIGNAL(textEdited(QString)), this, SLOT(setxtilt_offset_fromLineEdit())); 
    connect(this->ui->lineEdit_ytilt_offset, SIGNAL(textEdited(QString)), this, SLOT(setxtilt_offset_fromLineEdit())); 
    connect(this->ui->lineEdit_probeStepX, SIGNAL(textEdited(QString)), this, SLOT(setprobeStepX_fromLineEdit()));
    connect(this->ui->lineEdit_probeStepY, SIGNAL(textEdited(QString)), this, SLOT(setprobeStepY_fromLineEdit()));
    connect(this->ui->lineEdit_4Damax, SIGNAL(textEdited(QString)), this, SLOT(set4Damax_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowXMin, SIGNAL(textEdited(QString)), this, SLOT(setscan_WindowXMin_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowXMax, SIGNAL(textEdited(QString)), this, SLOT(setscan_WindowXMax_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowYMin, SIGNAL(textEdited(QString)), this, SLOT(setscan_WindowYMin_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowYMax, SIGNAL(textEdited(QString)), this, SLOT(setscan_WindowYMax_fromLineEdit()));
    connect(this->ui->lineEdit_scanWindowYMin, SIGNAL(editingFinished()), this, SLOT(setscan_WindowYMin_edited()));
    connect(this->ui->lineEdit_scanWindowYMax, SIGNAL(editingFinished()), this, SLOT(setscan_WindowYMax_edited()));
    connect(this->ui->lineEdit_scanWindowXMin, SIGNAL(textEdited(QString)), this, SLOT(updatePotentialDisplay()));
    connect(this->ui->lineEdit_scanWindowXMax, SIGNAL(textEdited(QString)), this, SLOT(updatePotentialDisplay()));
    connect(this->ui->lineEdit_scanWindowYMin, SIGNAL(textEdited(QString)), this, SLOT(updatePotentialDisplay()));
    connect(this->ui->lineEdit_scanWindowYMax, SIGNAL(textEdited(QString)), this, SLOT(updatePotentialDisplay()));
    connect(this->ui->lineEdit_interpFactor_y, SIGNAL(editingFinished()), this, SLOT(setinterpYSet_edited()));
    connect(this->ui->lineEdit_pixelSizeY, SIGNAL(editingFinished()), this, SLOT(setpixelSizeYSet_edited()));
    connect(this->ui->lineEdit_probeStepY, SIGNAL(editingFinished()), this, SLOT(setprobeStepYSet_edited()));
    connect(this->ui->lineEdit_probeTiltY, SIGNAL(editingFinished()), this, SLOT(setprobeTiltYSet_edited()));
    connect(this->ui->lineEdit_scanWindowXMin, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_scanWindowXMin()));
    connect(this->ui->lineEdit_scanWindowXMax, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_scanWindowXMax()));
    connect(this->ui->lineEdit_scanWindowYMin, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_scanWindowYMin()));
    connect(this->ui->lineEdit_scanWindowYMax, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_scanWindowYMax()));
    connect(this->ui->lineEdit_cellDimX, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_cellDimX()));
    connect(this->ui->lineEdit_cellDimY, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_cellDimY()));
    connect(this->ui->lineEdit_cellDimZ, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_cellDimZ()));
    connect(this->ui->lineEdit_tileX, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_tileX()));
    connect(this->ui->lineEdit_tileY, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_tileY()));
    connect(this->ui->lineEdit_tileZ, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_tileZ()));
    connect(this->ui->lineEdit_cellDimX, SIGNAL(textEdited(QString)), this, SLOT(userHasSetCellDims()));
    connect(this->ui->lineEdit_cellDimY, SIGNAL(textEdited(QString)), this, SLOT(userHasSetCellDims()));
    connect(this->ui->lineEdit_cellDimZ, SIGNAL(textEdited(QString)), this, SLOT(userHasSetCellDims()));
    connect(this->ui->lineEdit_pixelSizeX, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_pixelSizeX()));
    connect(this->ui->lineEdit_pixelSizeY, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_pixelSizeY()));
    connect(this->ui->lineEdit_interpFactor_x, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_interpFactor_x()));
    connect(this->ui->lineEdit_interpFactor_y, SIGNAL(editingFinished()), this, SLOT(checkInput_lineEdit_interpFactor_y()));
    connect(this->ui->lineEdit_E0, SIGNAL(textEdited(QString)), this, SLOT(setE0_fromLineEdit()));
    connect(this->ui->radBtn_PRISM, SIGNAL(clicked(bool)), this, SLOT(setAlgo_PRISM()));
    connect(this->ui->radBtn_Multislice, SIGNAL(clicked(bool)), this, SLOT(setAlgo_Multislice()));
    connect(this->ui->btn_calcPotential, SIGNAL(clicked(bool)), this, SLOT(calculatePotential()));
    connect(this->ui->btn_go, SIGNAL(clicked(bool)), this, SLOT(calculateAll()));
    connect(this->ui->btn_go_hrtem, SIGNAL(clicked(bool)), this, SLOT(calculateAllHRTEM()));
    connect(this->ui->lineEdit_slicemin, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits()));
    connect(this->ui->lineEdit_slicemax, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits()));
    connect(this->ui->slider_bothSlices, SIGNAL(valueChanged(int)), this, SLOT(moveBothPotentialSliders(int)));
    connect(this->ui->slider_bothDetectors, SIGNAL(valueChanged(int)), this, SLOT(moveBothDetectorSliders(int)));
    //    connect(this->ui->slider_slicemin,                 SIGNAL(valueChanged(int)),        this, SLOT(updateSlider_PotentialCombo(int)));
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
    connect(this->ui->lineEdit_contrast_outputMin, SIGNAL(editingFinished()), this, SLOT(updateContrastAngMin()));
    connect(this->ui->lineEdit_contrast_outputMax, SIGNAL(editingFinished()), this, SLOT(updateContrastAngMax()));
    connect(this->ui->lineEdit_contrastPotMin, SIGNAL(editingFinished()), this, SLOT(updateContrastPotMin()));
    connect(this->ui->lineEdit_contrastPotMax, SIGNAL(editingFinished()), this, SLOT(updateContrastPotMax()));
    connect(this->ui->tabWidget_2, SIGNAL(currentChanged(int)), this, SLOT(redrawImages()));
    connect(this->ui->tabWidget_3, SIGNAL(currentChanged(int)), this, SLOT(redrawImages()));
    connect(this->ui->btn_saveOutputImage, SIGNAL(clicked(bool)), this, SLOT(saveCurrentOutputImage()));
    connect(this->ui->comboBox_streamMode, SIGNAL(currentIndexChanged(int)), this, SLOT(setStreamingMode(int)));
    //connect(this->ui->checkBox_saveProjectedPotential, SIGNAL(toggled(bool)),            this, SLOT(toggleSaveProjectedPotential()));
    connect(this->ui->actionReset_Prismatic, SIGNAL(triggered()), this, SLOT(resetCalculation()));
    connect(this->ui->btn_calculateProbe, SIGNAL(clicked()), this, SLOT(calculateProbe()));
    connect(this->ui->actionReset_Prismatic, SIGNAL(triggered()), this, SLOT(resetLinks()));
    connect(this->ui->checkBox_2D, SIGNAL(toggled(bool)), this, SLOT(toggle2DOutput()));
    connect(this->ui->checkBox_3D, SIGNAL(toggled(bool)), this, SLOT(toggle3DOutput()));
    connect(this->ui->checkBox_4D, SIGNAL(toggled(bool)), this, SLOT(toggle4DOutput()));
    connect(this->ui->checkBox_crop4D, SIGNAL(toggled(bool)), this, SLOT(toggle4Dcrop()));
    connect(this->ui->checkBox_DPC_CoM, SIGNAL(toggled(bool)), this, SLOT(toggleDPC_CoM()));
    connect(this->ui->checkBox_PS, SIGNAL(toggled(bool)), this, SLOT(togglePotentialSlices()));
    connect(this->ui->checkBox_saveSMatrix, SIGNAL(toggled(bool)), this, SLOT(toggleSMatrixoutput()));
    connect(this->ui->checkBox_saveComplex, SIGNAL(toggled(bool)), this, SLOT(toggleComplexoutput()));
    connect(this->ui->checkBox_saveProbe, SIGNAL(toggled(bool)), this, SLOT(toggleProbeOutput()));
    connect(this->ui->checkBox_thermalEffects, SIGNAL(toggled(bool)), this, SLOT(toggleThermalEffects()));
    connect(this->ui->checkBox_matrixRefocus, SIGNAL(toggled(bool)), this, SLOT(togglematrixRefocus()));
    connect(this->ui->checkBox_potential3D, SIGNAL(toggled(bool)), this, SLOT(togglePotential3D()));
   // connect(this->ui->checkBox_occupancy, SIGNAL(toggled(bool)), this, SLOT(toggleOccupancy()));
    connect(this->ui->checkBox_NQS, SIGNAL(toggled(bool)), this, SLOT(toggleNyquist()));
    connect(this->ui->checkBox_sqrtIntensityPot, SIGNAL(toggled(bool)), this, SLOT(updatePotentialFloatImage()));
    connect(this->ui->checkBox_sqrtIntensityPot_2, SIGNAL(toggled(bool)), this, SLOT(updatePotentialFloatImage()));
    connect(this->ui->checkBox_log, SIGNAL(toggled(bool)), this, SLOT(updateProbeImages()));
    connect(this->ui->comboBox_colormap, SIGNAL(currentTextChanged(QString)), this, SLOT(changeColormap(QString)));
    connect(this->ui->comboBox_colormap,               SIGNAL(currentTextChanged(QString)), this, SLOT(updatePotentialFloatImage()));
    //    connect(this->ui->comboBox_colormap,               SIGNAL(currentTextChanged(QString)), this, SLOT(updateOutputFloatImage()));
    connect(this->ui->comboBox_colormap, SIGNAL(currentTextChanged(QString)), this, SLOT(updateAllImages()));

    //connections for the collapsable buttons
    connect(this->ui->btn_closebox, SIGNAL(clicked()), this, SLOT(collapseSample()));
    connect(this->ui->btn_closebox_2, SIGNAL(clicked()), this, SLOT(collapseSimulation()));
    connect(this->ui->btn_closebox_3, SIGNAL(clicked()), this, SLOT(collapseStem()));
    connect(this->ui->btn_closebox_4, SIGNAL(clicked()), this, SLOT(collapseHrtem()));
    connect(this->ui->btn_closebox_5, SIGNAL(clicked()), this, SLOT(collapseOutput()));
    connect(this->ui->btn_closebox_7, SIGNAL(clicked()), this, SLOT(collapseComputational()));


    //connections for changing the theme/field of the application
    connect(this->ui->actionDarkField, SIGNAL(triggered()), this, SLOT(darkField()));
    connect(this->ui->actionLightField, SIGNAL(triggered()), this, SLOT(lightField()));







    //    connect(this->ui->tabs,                            SIGNAL(currentChanged(int)),this, SLOT(updatePotentialDisplay()));
    updateAlphaMax();
    //    ui->lbl_image_potential->setPixmap(QPixmap::fromImage(potentialImage.scaled(ui->tabs->width(),
    //                                                                                ui->tabs->height(),
    //                                                                                Qt::KeepAspectRatio)));
}



void PRISMMainWindow::updateDisplay(){
    std::stringstream ss;
	ss << this->meta->realspacePixelSize[0];
	this->ui->lineEdit_pixelSizeY->setText(QString::fromStdString(ss.str()));
	ss.str("");
	ss << this->meta->realspacePixelSize[1];
	this->ui->lineEdit_pixelSizeX->setText(QString::fromStdString(ss.str()));
	ss.str("");
    ss << this->meta->interpolationFactorX;
    this->ui->lineEdit_interpFactor_x->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << this->meta->interpolationFactorY;
    this->ui->lineEdit_interpFactor_y->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << this->meta->potBound;
    this->ui->lineEdit_potbound->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeSemiangle * 1e3);
    this->ui->lineEdit_probeSemiangle->setText(QString::fromStdString(ss.str()));
    ss.str("");
	ss << (this->meta->zStart);
    this->ui->lineEdit_zStart->setText(QString::fromStdString(ss.str()));
    ss.str("");
    // ss << (this->meta->alphaBeamMax * 1e3);
    // this->ui->lineEdit_alphaBeamMax->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->integrationAngleMin * 1e3);
    this->ui->lineEdit_2D_inner->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->integrationAngleMax * 1e3);
    this->ui->lineEdit_2D_outer->setText(QString::fromStdString(ss.str()));
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
    ss << (this->meta->probeStepX);
    this->ui->lineEdit_probeStepX->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeStepY);
    this->ui->lineEdit_probeStepY->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeDefocus);
    this->ui->lineEdit_probeDefocus->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeDefocus_min);
    this->ui->lineEdit_dfr_min->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeDefocus_max);
    this->ui->lineEdit_dfr_max->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeDefocus_step);
    this->ui->lineEdit_dfr_step->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeDefocus_sigma);
    this->ui->lineEdit_dfs->setText(QString::fromStdString(ss.str()));
    ss.str("");
	ss << (this->meta->C3);
	this->ui->lineEdit_C3->setText(QString::fromStdString(ss.str()));
	ss.str("");
	ss << (this->meta->C5);
	this->ui->lineEdit_C5->setText(QString::fromStdString(ss.str()));
	ss.str("");
    ss << (this->meta->probeXtilt * 1e3);
    this->ui->lineEdit_probeTiltX->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->probeYtilt * 1e3);
    this->ui->lineEdit_probeTiltY->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->minXtilt * 1e3);
    this->ui->lineEdit_xtt_min->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->maxXtilt * 1e3);
    this->ui->lineEdit_xtt_max->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->xTiltStep * 1e3);
    this->ui->lineEdit_xtt_step->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->minYtilt * 1e3);
    this->ui->lineEdit_ytt_min->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->maxYtilt * 1e3);
    this->ui->lineEdit_ytt_max->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->yTiltStep * 1e3);
    this->ui->lineEdit_ytt_step->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->minRtilt * 1e3);
    this->ui->lineEdit_rtt_min->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->maxRtilt * 1e3);
    this->ui->lineEdit_rtt_max->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->xTiltOffset * 1e3);
    this->ui->lineEdit_xtilt_offset->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->yTiltOffset * 1e3);
    this->ui->lineEdit_ytilt_offset->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->detectorAngleStep * 1e3);
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
    ss << (this->meta->randomSeed);
    this->ui->lineEdit_randomSeed->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->crop4Damax);
    this->ui->lineEdit_4Damax->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->batchSizeTargetCPU);
    this->ui->lineEdit_batchCPU->setText(QString::fromStdString(ss.str()));
    ss.str("");
    ss << (this->meta->batchSizeTargetGPU);
    this->ui->lineEdit_batchGPU->setText(QString::fromStdString(ss.str()));
    ss.str("");

    this->ui->lineEdit_scanWindowXMin->setCursorPosition(0);
    this->ui->lineEdit_scanWindowXMax->setCursorPosition(0);
    this->ui->lineEdit_scanWindowYMin->setCursorPosition(0);
    this->ui->lineEdit_scanWindowYMax->setCursorPosition(0);
    this->ui->lineEdit_randomSeed->setCursorPosition(0);
    ui->lineEdit_4Damax->setEnabled(false);

    this->ui->lineEdit_outputfile->setText(QString::fromStdString(ss.str()));
    this->ui->spinBox_numGPUs->setValue(this->meta->numGPUs);
    this->ui->spinBox_numThreads->setValue(this->meta->numThreads);
    this->ui->spinBox_numFP->setValue(this->meta->numFP);
    this->ui->spinBox_numNS->setValue(this->meta->numSlices);
    this->ui->spinBox_zSampling->setValue(this->meta->zSampling);
    this->ui->spinBox_numStreams->setValue(this->meta->numStreamsPerGPU);
    ui->checkBox_potential3D->setChecked(meta->potential3D);
    ui->checkBox_matrixRefocus->setChecked(meta->matrixRefocus);
    ui->checkBox_thermalEffects->setChecked(meta->includeThermalEffects);
   // ui->checkBox_occupancy->setChecked(meta->includeOccupancy);
    ui->checkBox_NQS->setChecked(meta->nyquistSampling);
    ui->checkBox_2D->setChecked(meta->save2DOutput);
    ui->checkBox_3D->setChecked(meta->save3DOutput);
    ui->checkBox_4D->setChecked(meta->save4DOutput);
    ui->checkBox_DPC_CoM->setChecked(meta->saveDPC_CoM);
    ui->checkBox_PS->setChecked(meta->savePotentialSlices);
    ui->checkBox_saveSMatrix->setChecked(meta->saveSMatrix);
    ui->checkBox_saveComplex->setChecked(meta->saveComplexOutputWave);
    ui->checkBox_saveProbe->setChecked(meta->saveProbe);


    switch (this->meta->algorithm){
        case Prismatic::Algorithm::PRISM :
            this->ui->radBtn_PRISM->setChecked(true);
            this->ui->radBtn_Multislice->setChecked(false);
            break;
        case Prismatic::Algorithm::Multislice : 
             this->ui->radBtn_PRISM->setChecked(false);
            this->ui->radBtn_Multislice->setChecked(true);
            break;
    }
#ifndef PRISMATIC_ENABLE_GPU
    this->ui->spinBox_numGPUs->setEnabled(false);
    this->ui->spinBox_numStreams->setEnabled(false);
    this->ui->lineEdit_batchGPU ->setEnabled(false);
    this->ui->comboBox_streamMode->setEnabled(false);
#endif //PRISMATIC_ENABLE_GPU


    ui->lbl_angstrom->setText(QString::fromUtf8("\u212B"));
    ui->lbl_sliceThickness->setText(QString::fromUtf8("Slice\nThickness (\u212B)"));
    ui->lbl_probeStep->setText(QString::fromUtf8("Probe Step (\u212B)"));
    ui->lbl_alphaMax->setText(QString::fromUtf8("\u03B1 max = ??"));
    // ui->lbl_alphaBeamMax->setText(QString::fromUtf8("Probe \u03B1 limit (mrads)"));
    ui->lbl_lambda->setText(QString::fromUtf8("\u03BB = ") + QString::number(calculateLambda(*meta)) + QString::fromUtf8("\u212B"));
    ui->lbl_potBound->setText(QString::fromUtf8("Potential\nBound (\u212B)"));
    ui->lbl_pixelSize->setText(QString::fromUtf8("Pixel\nSize (\u212B)"));
    ui->lbl_defocus->setText(QString::fromUtf8("C1 (defocus)(\u212B)"));
    ui->lbl_C3->setText(QString::fromUtf8("C3 (\u212B)"));
    ui->lbl_C5->setText(QString::fromUtf8("C5 (\u212B)"));
    ui->label_Xprobe->setText(QString::fromUtf8("X (\u212B)"));
    ui->label_Yprobe->setText(QString::fromUtf8("Y (\u212B)"));

    this->ui->lineEdit_outputfile->setText(QString::fromStdString(this->meta->filenameOutput));
}


void PRISMMainWindow::readParams(std::string param_filename){
    if(!Prismatic::parseParamFile(*this->meta, param_filename)){
        displayErrorReadingParamsDialog();
    } else {
        if (validateFilename(this->meta->filenameAtoms)){
            updateUCdims(this->meta->filenameAtoms);
        }
    }
    updateDisplay();
}

void PRISMMainWindow::selectParameterFile(){
    QString filename;
    filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename, tr("Parameter File(*.txt);;All files(*)"));
    if (validateFilename(filename.toStdString())){
        readParams(filename.toStdString());
    }
}

void PRISMMainWindow::writeParameterFile(){
    QString filename;
    filename = QFileDialog::getSaveFileName(this, tr("ExistingFile"), filename, tr("Parameter File(*.txt);;All files(*)"));
    if (validateWriteFilename(filename.toStdString())){
        Prismatic::writeParamFile(*this->meta, filename.toStdString());
    }
}

void PRISMMainWindow::readAberrationFile(){
    QString filename;
    filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename, tr("Aberration File(*.txt);;All files(*)"));
    if (validateFilename(filename.toStdString())){
        this->meta->aberrations = Prismatic::readAberrations(filename.toStdString());
        this->meta->arbitraryAberrations = true;
    }
}

void PRISMMainWindow::readProbeFile(){
    QString filename;
    filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename, tr("Probe File(*.txt);;All files(*)"));
    if (validateFilename(filename.toStdString())){
        Prismatic::readProbes(filename.toStdString(), this->meta->probes_x, this->meta->probes_y);
        this->meta->arbitraryProbes = true;
    }
}

void PRISMMainWindow::setAlgo_PRISM(){
	std::cout << "Setting algorithm to PRISM" << std::endl;
	setAlgo(Prismatic::Algorithm::PRISM);
    ui->lineEdit_interpFactor_x->setEnabled(true);
    ui->lineEdit_interpFactor_y->setEnabled(true);
    // ui->lineEdit_alphaBeamMax->setEnabled(true);
}

void PRISMMainWindow::setAlgo_Multislice(){
	std::cout << "Setting algorithm to Multislice" << std::endl;
	setAlgo(Prismatic::Algorithm::Multislice);
    ui->lineEdit_interpFactor_x->setDisabled(true);
    ui->lineEdit_interpFactor_y->setDisabled(true);
    // ui->lineEdit_alphaBeamMax->setDisabled(true);
}

void PRISMMainWindow::setAlgo(const Prismatic::Algorithm algo){
	this->meta->algorithm = algo;
    resetCalculation();
}

void PRISMMainWindow::setInterpolationFactorX(){
	bool flag;
    const size_t& new_f = this->ui->lineEdit_interpFactor_x->text().toUInt(&flag);
	if (flag){
		std::cout << "Setting interpolation factor X to " << new_f << std::endl;
		this->meta->interpolationFactorX = new_f;
        if (!interpYSet){
            ui->lineEdit_interpFactor_y->setText(QString::number(new_f));
            std::cout << "Setting interpolation factor Y to " << new_f << std::endl;
            this->meta->interpolationFactorY = new_f;
//            setInterpolationFactorY();
        }
    }
    resetCalculation();
}

void PRISMMainWindow::setInterpolationFactorY(){
    bool flag;
    const size_t& new_f = this->ui->lineEdit_interpFactor_y->text().toUInt(&flag);
    if (flag){
        std::cout << "Setting interpolation factor Y to " << new_f << std::endl;
        this->meta->interpolationFactorY = new_f;
    }
    resetCalculation();
}

void PRISMMainWindow::setFilenameAtoms_fromDialog(){
	QString filename;
    filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename, tr("Atomic Model(*.xyz *.XYZ);;All files(*)"));
    if (validateFilename(filename.toStdString())){
        updateUCdims(filename.toStdString());
        meta->userSpecifiedCelldims = false;
    }
    resetCalculation();
}

void PRISMMainWindow::updateUCdims(const std::string& filename){
    // get the unit cell dimensions from the input file (if possible)
    bool error_reading = false;
    std::array<double, 3> uc_dims;
    try {
        uc_dims = Prismatic::peekDims_xyz(filename);
    } catch (...){
        error_reading = true;
    }
    if (error_reading){
        displayErrorReadingAtomsDialog();
    }else{
        this->setFilenameAtoms(filename);
        ui->btn_go->setEnabled(true);
        ui->btn_go_hrtem->setEnabled(true);
        ui->btn_calcPotential->setEnabled(true);
        this->setWindowTitle(QString::fromStdString(std::string("Prismatic (") + std::string(filename + std::string(")"))));
        if (uc_dims[0]>0){
            // update gui
            ui->lineEdit_cellDimX->setText(QString::number(uc_dims[2]));
            ui->lineEdit_cellDimY->setText(QString::number(uc_dims[1]));
            ui->lineEdit_cellDimZ->setText(QString::number(uc_dims[0]));

            // move cursor of the cell dimension line edits back to 0 so they are easy to read
            ui->lineEdit_cellDimX->setCursorPosition(0);
            ui->lineEdit_cellDimY->setCursorPosition(0);
            ui->lineEdit_cellDimZ->setCursorPosition(0);

            meta->cellDim[0] = uc_dims[0];
            meta->cellDim[1] = uc_dims[1];
            meta->cellDim[2] = uc_dims[2];
        }
    }
}

void PRISMMainWindow::setFilenameOutput_fromDialog(){
	QString filename;
	filename = QFileDialog::getOpenFileName(this, tr("AnyFile"), filename);
    this->ui->lineEdit_outputfile->setText(filename);
	this->setFilenameOutput(filename.toStdString());
    resetCalculation();
}


void PRISMMainWindow::setFilenameOutput_fromLineEdit(){
    const std::string& filename = this->ui->lineEdit_outputfile->text().toStdString();
	this->setFilenameOutput(filename);
    resetCalculation();
}

void PRISMMainWindow::setFilenameAtoms(const std::string& filename){
	std::cout << "Setting atoms filename to " << filename << std::endl;
	this->meta->filenameAtoms = filename;
    resetCalculation();
}

void PRISMMainWindow::setFilenameOutput(const std::string& filename){
	std::cout << "Setting output filename to " << filename << std::endl;
	this->meta->filenameOutput = filename;
    resetCalculation();
}


void PRISMMainWindow::setNumGPUs(const int& num){
    if (num >= 0){
        this->meta->numGPUs = num;
        std::cout << "Setting number of GPUs to " << num << std::endl;
        QMutexLocker gatekeeper(&dataLock);
        this->pars.meta.numGPUs = num;
    }
}

void PRISMMainWindow::setNumThreads(const int& num){
    std::cout << "Also do CPU work " << pars.meta.alsoDoCPUWork << std::endl;
    if (num > 0){
        this->meta->numThreads = num;
        std::cout << "Setting number of CPU Threads to " << num << std::endl;
        QMutexLocker gatekeeper(&dataLock);
        this->meta->alsoDoCPUWork=true;
        this->pars.meta.numThreads = num;
    }
    else if(num==0){
      this->meta->numThreads = 1;
      std::cout << "Setting number of CPU threads to " << num <<" (GPU only calculation)"<< std::endl;
      this->meta->alsoDoCPUWork=false;
      this->pars.meta.numThreads = 1;
    }
    std::cout << "Also do CPU work " << pars.meta.alsoDoCPUWork << std::endl;
}

void PRISMMainWindow::setNumStreams(const int& num){
    if (num >= 0){
        this->meta->numStreamsPerGPU = num;
        std::cout << "Setting number of CUDA streams per GPU to " << num << std::endl;
        QMutexLocker gatekeeper(&dataLock);
        this->pars.meta.numStreamsPerGPU = num;
    }
}

void PRISMMainWindow::setNumFP(const int& num){
    if (num > 0){
        this->meta->numFP = num;
        std::cout << "Setting number of frozen phonon configurations to " << num << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setNumNS(const int& num){
    if (num > 0){
        this->meta->numSlices = num;
        std::cout << "Setting number of slices for intermediate output steps to " << num << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setzSampling(const int& num){
    if (num > 0){
        this->meta->zSampling = num;
        std::cout << "Setting number of subslices for 3D potential integration to " << num << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setPixelSizeX_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val =(PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_pixelSizeX->text().toDouble(&flag);
    if (flag & (val > 0)){
        this->meta->realspacePixelSize[1] = val;
        std::cout << "Setting X pixel size to " << val << " Angstroms" << std::endl;
        if (!pixelSizeYSet){
            ui->lineEdit_pixelSizeY->setText(QString::number(val));
            this->meta->realspacePixelSize[0] = val;
            std::cout << "Setting Y pixel size to " << val << " Angstroms" << std::endl;
//            setPixelSizeY_fromLineEdit();
        }
        updateAlphaMax();
    }
    resetCalculation();
}

void PRISMMainWindow::setPixelSizeY_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val =(PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_pixelSizeY->text().toDouble(&flag);
    if (flag & (val > 0)){
        this->meta->realspacePixelSize[0] = val;
        std::cout << "Setting Y pixel size to " << val << " Angstroms" << std::endl;
        updateAlphaMax();
    }
    resetCalculation();
}


void PRISMMainWindow::setPotBound_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_potbound->text().toDouble(&flag);
    if (flag){
        this->meta->potBound = val;
        std::cout << "Setting potential bound to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobeSemiangle_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_probeSemiangle->text().toDouble(&flag);
    if (flag){
        this->meta->probeSemiangle = val / 1000;
        std::cout << "Setting probe semiangle to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setzStart_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_zStart->text().toDouble(&flag);
    if (flag){
        this->meta->zStart = val;
        std::cout << "Setting intermediate output to begin after" << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

// void PRISMMainWindow::setalphaBeamMax_fromLineEdit(){
//     bool flag = false;
//     PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_alphaBeamMax->text().toDouble(&flag);
//     if (flag){
//         this->meta->alphaBeamMax = val / 1000;
//         std::cout << "Setting maximum PRISM probe scattering angle to " << val << " mrad" << std::endl;
//     }
//     resetCalculation();
// }

void PRISMMainWindow::set2D_innerAngle_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_2D_inner->text().toDouble(&flag);
    if (flag){
        this->meta->integrationAngleMin = val /1000;
        std::cout << "Setting annular detector inner angle to" << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::set2D_outerAngle_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_2D_inner->text().toDouble(&flag);
    if (flag){
        this->meta->integrationAngleMax = val /1000;
        std::cout << "Setting annular detector outer angle to" << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setSliceThickness_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_sliceThickness->text().toDouble(&flag);
    if (flag){
        this->meta->sliceThickness = val;
        std::cout << "Setting sliceThickness to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setCellDimX_fromLineEdit(){
    bool flag = false;
    double val = this->ui->lineEdit_cellDimX->text().toDouble(&flag);
    if (flag){
        this->meta->cellDim[2] = (PRISMATIC_FLOAT_PRECISION)val;
        std::cout << "Setting X cell dimension to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setCellDimY_fromLineEdit(){
    bool flag = false;
    double val = this->ui->lineEdit_cellDimY->text().toDouble(&flag);
    if (flag){
        this->meta->cellDim[1] = (PRISMATIC_FLOAT_PRECISION)val;
        std::cout << "Setting Y cell dimension to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setCellDimZ_fromLineEdit(){
    bool flag = false;
    double val = this->ui->lineEdit_cellDimZ->text().toDouble(&flag);
    if (flag){
        this->meta->cellDim[0] = (PRISMATIC_FLOAT_PRECISION)val;
        std::cout << "Setting Z cell dimension to " << val << " Angstroms" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setTileX_fromLineEdit(){
    bool flag = false;
    int val = this->ui->lineEdit_tileX->text().toInt(&flag);
    if (flag){
        this->meta->tileX = (size_t)val;
        std::cout << "Setting tileX to " << val << " UCs" << std::endl;
        updateAlphaMax();
    }
    resetCalculation();
}

void PRISMMainWindow::setTileY_fromLineEdit(){
    bool flag = false;
    int val = this->ui->lineEdit_tileY->text().toInt(&flag);
    if (flag){
        this->meta->tileY = (size_t)val;
        std::cout << "Setting tileY to " << val << " UCs" << std::endl;
        updateAlphaMax();
    }
    resetCalculation();
}

void PRISMMainWindow::setTileZ_fromLineEdit(){
    bool flag = false;
    int val = this->ui->lineEdit_tileZ->text().toInt(&flag);
    if (flag){
        this->meta->tileZ = (size_t)val;
        std::cout << "Setting tileZ to " << val << " UCs" << std::endl;
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setRandomSeed_fromLineEdit(){
    bool flag = false;
    int val = this->ui->lineEdit_randomSeed->text().toInt(&flag);
    if (flag){
        this->meta->randomSeed = (size_t)val;
        std::cout << "Setting random seed to " << val << std::endl;
    }
    resetCalculation();
}


void PRISMMainWindow::setE0_fromLineEdit(){
    bool flag = false;
    double val = this->ui->lineEdit_E0->text().toDouble(&flag);
    if (flag){
        meta->E0 = val * 1e3;
        std::cout << "Setting E0 to " << val << " keV" << std::endl;
        ui->lbl_lambda->setText(QString::fromUtf8("\u03BB = ") + QString::number(calculateLambda(*meta)) + QString::fromUtf8("\u212B"));
        updateAlphaMax();

    }
    resetCalculation();
}

void PRISMMainWindow::setprobeStepX_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_probeStepX->text().toDouble(&flag);
    if (flag & (val > 0)){
        this->meta->probeStepX = val;
        std::cout << "Setting probeStepX to " << val << " Angstroms" << std::endl;
        if (!probeStepYSet){
            this->ui->lineEdit_probeStepY->setText(this->ui->lineEdit_probeStepX->text());
            this->meta->probeStepY = val;
            std::cout << "Setting probeStepY to " << val << " Angstroms" << std::endl;
        }
    }
    resetCalculation();
}

void PRISMMainWindow::setprobeStepY_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_probeStepY->text().toDouble(&flag);
    if (flag & (val > 0)){
        this->meta->probeStepY = val;
        std::cout << "Setting probeStepY to " << val << " Angstroms" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_defocus_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_probeDefocus->text().toDouble(&flag);
    if (flag){
        this->meta->probeDefocus = val;
        std::cout << "Setting probe defocus to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::set_dfr_min_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_dfr_min->text().toDouble(&flag);
    if (flag){
        this->meta->probeDefocus_min = val;
        std::cout << "Setting probe defocus range min to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::set_dfr_max_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_dfr_max->text().toDouble(&flag);
    if (flag){
        this->meta->probeDefocus_max = val;
        std::cout << "Setting probe defocus range max to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::set_dfr_step_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_dfr_step->text().toDouble(&flag);
    if (flag){
        this->meta->probeDefocus_step = val;
        std::cout << "Setting probe defocus range step to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::set_dfs_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_dfs->text().toDouble(&flag);
    if (flag){
        this->meta->probeDefocus_sigma = val;
        std::cout << "Setting probe defocus range step to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_C3_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_C3->text().toDouble(&flag);
    if (flag){
        this->meta->C3 = val;
        std::cout << "Setting C3 to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_C5_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_C5->text().toDouble(&flag);
    if (flag){
        this->meta->C5 = val;
        std::cout << "Setting C5 to " << val << " Angstroms" <<  std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setdetectorAngleStep_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_detectorAngle->text().toDouble(&flag);
    if (flag){
        this->meta->detectorAngleStep = val / 1000;
        std::cout << "Setting detector angle step to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}


void PRISMMainWindow::setprobe_Xtilt_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_probeTiltX->text().toDouble(&flag);
    if (flag){
        this->meta->probeXtilt = val / 1000;
        std::cout << "Setting probe X tilt to " << val << " mrad" << std::endl;
        if (!probeTiltYSet){
            ui->lineEdit_probeTiltY->setText(ui->lineEdit_probeTiltX->text());
            this->meta->probeYtilt = val / 1000;
            std::cout << "Setting probe Y tilt to " << val << " mrad" << std::endl;
        }
    }
    resetCalculation();
}

void PRISMMainWindow::setprobe_Ytilt_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_probeTiltY->text().toDouble(&flag);
    if (flag){
        this->meta->probeYtilt =  val / 1000;
        std::cout << "Setting probe Y tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setxtt_min_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_xtt_min->text().toDouble(&flag);
    if (flag){
        this->meta->minXtilt =  val / 1000;
        std::cout << "Setting HRTEM min X tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setxtt_max_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_xtt_max->text().toDouble(&flag);
    if (flag){
        this->meta->maxXtilt =  val / 1000;
        std::cout << "Setting HRTEM max X tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setxtt_step_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_xtt_step->text().toDouble(&flag);
    if (flag){
        this->meta->xTiltStep =  val / 1000;
        std::cout << "Setting HRTEM X tilt step to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setytt_min_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_ytt_min->text().toDouble(&flag);
    if (flag){
        this->meta->minXtilt =  val / 1000;
        std::cout << "Setting HRTEM min Y tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setytt_max_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_ytt_max->text().toDouble(&flag);
    if (flag){
        this->meta->maxXtilt =  val / 1000;
        std::cout << "Setting HRTEM max Y tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setytt_step_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_ytt_step->text().toDouble(&flag);
    if (flag){
        this->meta->yTiltStep =  val / 1000;
        std::cout << "Setting HRTEM Y tilt step to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setrtt_min_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_rtt_min->text().toDouble(&flag);
    if (flag){
        this->meta->minXtilt =  val / 1000;
        std::cout << "Setting HRTEM min R tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setrtt_max_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_rtt_max->text().toDouble(&flag);
    if (flag){
        this->meta->maxXtilt =  val / 1000;
        std::cout << "Setting HRTEM max R tilt to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setxtilt_offset_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_xtilt_offset->text().toDouble(&flag);
    if (flag){
        this->meta->xTiltOffset =  val / 1000;
        std::cout << "Setting HRTEM X tilt offset to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setytilt_offset_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_ytilt_offset->text().toDouble(&flag);
    if (flag){
        this->meta->yTiltOffset =  val / 1000;
        std::cout << "Setting HRTEM Y tilt offset to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setBatchGPU_fromLineEdit(){
    bool flag = false;
    int val = this->ui->lineEdit_batchGPU->text().toInt(&flag);
    if (flag){
        this->meta->batchSizeTargetGPU = val;
        this->meta->batchSizeGPU = val;
        QMutexLocker gatekeeper(&dataLock);
        this->pars.meta.batchSizeTargetGPU = val;
        this->pars.meta.batchSizeGPU = val;
        std::cout << "Setting batch size (GPU) to " << val << std::endl;
    }
}

void PRISMMainWindow::setBatchCPU_fromLineEdit(){
    bool flag = false;
    int val = this->ui->lineEdit_batchCPU->text().toInt(&flag);
    if (flag){
        this->meta->batchSizeTargetCPU = val;
        this->meta->batchSizeCPU = val;
        QMutexLocker gatekeeper(&dataLock);
        this->pars.meta.batchSizeTargetCPU = val;
        this->pars.meta.batchSizeCPU = val;
        std::cout << "Setting batch size (CPU) to " << val << std::endl;
    }
}

void PRISMMainWindow::setscan_WindowXMin_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_scanWindowXMin->text().toDouble(&flag);
    if (flag){
        val = std::min(this->meta->scanWindowXMax, val);
        this->meta->scanWindowXMin = val;
        std::cout << "Setting scan window X min to " << val << std::endl;
        if (!minWindowYSet){
            this->ui->lineEdit_scanWindowYMin->setText(QString::number(val));
            this->meta->scanWindowYMin = val;
            std::cout << "Setting scan window Y min to " << val << std::endl;
        }
    }
    resetCalculation();
}

void PRISMMainWindow::setscan_WindowXMax_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_scanWindowXMax->text().toDouble(&flag);
    if (flag){
        val = std::max(this->meta->scanWindowXMin, val);
        this->meta->scanWindowXMax = val;
        std::cout << "Setting scan window X max to " << val << std::endl;
        if (!maxWindowYSet){
            this->ui->lineEdit_scanWindowYMax->setText(QString::number(val));
            this->meta->scanWindowYMax = val;
            std::cout << "Setting scan window Y max to " << val << std::endl;
        }
    }
    resetCalculation();
}

void PRISMMainWindow::setscan_WindowYMin_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_scanWindowYMin->text().toDouble(&flag);
    if (flag){
        val = std::min(this->meta->scanWindowYMax, val);
        this->meta->scanWindowYMin = val;
        std::cout << "Setting scan window Y min to " << val << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::setscan_WindowYMax_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_scanWindowYMax->text().toDouble(&flag);
    if (flag){
        val = std::max(this->meta->scanWindowYMin, val);
        this->meta->scanWindowYMax = val;
        std::cout << "Setting scan window Y max to " << val << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::calculatePotential(){
    resetCalculation();
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();
    PotentialThread *worker = new PotentialThread(this, progressbar);
    worker->meta.toString();
    connect(worker, SIGNAL(signalErrorReadingAtomsDialog()), this, SLOT(displayErrorReadingAtomsDialog()));
    //connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
    connect(worker, SIGNAL(finished()), this, SLOT(updatePotentialImage()));
    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
    worker->start();
}

//void PRISMMainWindow::calculateSMatrix(){
//    prism_progressbar *progressbar = new prism_progressbar(this);
//    progressbar->show();
//    SMatrixThread *worker = new SMatrixThread(this, progressbar);
//    worker->meta.toString();
////    connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
//    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
//    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
//    worker->start();
//}

void PRISMMainWindow::openSaveAtomsDialog(){
    SaveAtomicCoordinatesDialog *dialog = new SaveAtomicCoordinatesDialog(this);
    //dialog->setFilenameText(QString::fromStdString(meta->filenameAtoms.substr(0, meta->filenameAtoms.rfind("."))));
    std::stringstream tileX_ss; tileX_ss << meta->tileX;
    std::stringstream tileY_ss; tileY_ss << meta->tileY;
    std::stringstream tileZ_ss; tileZ_ss << meta->tileZ;
    dialog->setFilenameText(QString::fromStdString(
		    meta->filenameAtoms.substr(0, meta->filenameAtoms.rfind(".")) +
            std::string("tiled") +
            tileX_ss.str() +
            std::string("x") +
            tileY_ss.str()  +
            std::string("x") +
            tileZ_ss.str()  +
            std::string(".XYZ")));

    dialog->setCommentText(QString::fromStdString(
		    meta->filenameAtoms +
		    std::string(" tiled ") +
		    tileX_ss.str() +
		    std::string("x") +
		    tileY_ss.str()  +
		    std::string("x") +
		    tileZ_ss.str() ));

    dialog->show();
    connect(dialog, SIGNAL(accepted()), dialog, SLOT(SaveAtomCoords()));
    connect(dialog, SIGNAL(signalSaveAtomCoords(QString, QString)), this, SLOT(saveAtomCoords(QString, QString)));
    connect(dialog, SIGNAL(accepted()), dialog, SLOT(deleteLater()));
}

void PRISMMainWindow::saveAtomCoords(QString filename, QString comment){
//    std::cout << "Saving tiled coords from file " << meta->filenameAtoms << std::endl;
//    std::cout << "Saving tiled coords with filename " << filename.toStdString() << std::endl;
//    std::cout << "Saving tiled coords with comment " << comment.toStdString() << std::endl;
    bool error_reading = false;
    Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars;
    try {
        pars = Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION>(*meta);
    } catch(const std::runtime_error &e){
        error_reading = true;
    }
    if (error_reading){
        displayErrorReadingAtomsDialog();
    }else{
        Prismatic::to_xyz(pars.atoms, filename.toStdString(), comment.toStdString(), pars.tiledCellDim[2], pars.tiledCellDim[1], pars.tiledCellDim[0]);
//        Prismatic::to_xyz(pars.atoms, filename.toStdString(), comment.toStdString(), meta->cellDim[2], meta->cellDim[1], meta->cellDim[0]);

    }
}

void PRISMMainWindow::calculateAll(){
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();
    this->setFilenameOutput_fromLineEdit();

    FullCalcThread *worker = new FullCalcThread(this, progressbar);
    worker->meta.toString();
    connect(worker, SIGNAL(signalErrorReadingAtomsDialog()), this, SLOT(displayErrorReadingAtomsDialog()));
    connect(worker, SIGNAL(overwriteWarning()),this,SLOT(preventOverwrite()),Qt::BlockingQueuedConnection);
    connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));
    if(this->meta->algorithm == Prismatic::Algorithm::HRTEM){
        connect(worker, SIGNAL(outputCalculated()), this, SLOT(updateOutputImage_HRTEM()));
    }else{
        connect(worker, SIGNAL(outputCalculated()), this, SLOT(updateOutputImage()));
    }
    connect(worker, SIGNAL(outputCalculated()), this, SLOT(enableOutputWidgets()));
    connect(worker, SIGNAL(signalTitle(const QString)), progressbar, SLOT(setTitle(const QString)));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(close()));
    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
    worker->start();
    Prismatic::writeParamFile(*this->meta, get_default_parameter_filename());
}

void PRISMMainWindow::calculateAllHRTEM(){
    Prismatic::Algorithm stem_algo = this->meta->algorithm;

    this->meta->algorithm = Prismatic::Algorithm::HRTEM;
    calculateAll();

    setAlgo(stem_algo);
}

void PRISMMainWindow::preventOverwrite(){
    QMessageBox *msgBox = new QMessageBox(this);
    QString msgTxt = "Output file already exists. Do you want to overwrite the file?";
    msgBox->setText(msgTxt);
    msgBox->setStandardButtons(QMessageBox::No | QMessageBox::Yes);
    msgBox->show();
    int ret;
    ret = msgBox->exec();
    switch (ret) 
    {
        case QMessageBox::Yes:
            //delete file and move on
            //remove(params.meta.filenameOutput.c_str());
            flipOverwrite();
            break;
        case QMessageBox::No:
            //exit from thread
            break;
    }
}

void PRISMMainWindow::calculateProbe(){
    prism_progressbar *progressbar = new prism_progressbar(this);
    progressbar->show();
    bool flagX = false;
    bool flagY = false;
    PRISMATIC_FLOAT_PRECISION X = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_probeX->text().toDouble(&flagX);
    PRISMATIC_FLOAT_PRECISION Y = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_probeY->text().toDouble(&flagY);
    if (flagX & flagY){
        currently_calculated_X = X;
        currently_calculated_Y = Y;
        ProbeThread *worker = new ProbeThread(this, X, Y, progressbar, ui->checkBox_log->isChecked());
        connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
        connect(worker, SIGNAL(finished()), this, SLOT(updatePotentialImage()));
        connect(worker, SIGNAL(signalErrorReadingAtomsDialog()), this, SLOT(displayErrorReadingAtomsDialog()));
        connect(worker, SIGNAL(signal_RReal(QString)), this, SLOT(update_RReal(QString)));
        connect(worker, SIGNAL(signal_RK(QString)), this, SLOT(update_RK(QString)));
        connect(worker, SIGNAL(finished()), progressbar, SLOT(close()));
        connect(worker, SIGNAL(finished()), progressbar, SLOT(deleteLater()));
        connect(worker, SIGNAL(potentialCalculated()), this, SLOT(updatePotentialImage()));

        connect(worker, SIGNAL(signal_pearsonReal(QString)), this, SLOT(update_pearsonReal(QString)));
        connect(worker, SIGNAL(signal_pearsonK(QString)), this, SLOT(update_pearsonK(QString)));

//        connect(worker, SIGNAL(signalProbeK_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
//                this, SLOT(updateProbeK_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
//        connect(worker, SIGNAL(signalProbeR_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
//                this, SLOT(updateProbeR_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
//        connect(worker, SIGNAL(signalProbeK_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
//                this, SLOT(updateProbeK_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
//        connect(worker, SIGNAL(signalProbeR_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
//                this, SLOT(updateProbeR_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
        connect(worker, SIGNAL(signalProbeK_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
                this, SLOT(probeK_PRISMReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
        connect(worker, SIGNAL(signalProbeR_PRISM(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
                this, SLOT(probeR_PRISMReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
        connect(worker, SIGNAL(signalProbeK_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
                this, SLOT(probeK_MultisliceReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
        connect(worker, SIGNAL(signalProbeR_Multislice(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
                this, SLOT(probeR_MultisliceReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));


        connect(worker, SIGNAL(signalProbe_diffR(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
                this, SLOT(calculateProbe_diffR(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
        connect(worker, SIGNAL(signalProbe_diffK(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)),
                this, SLOT(calculateProbe_diffK(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION>)));
        worker->start();
    }
}

void PRISMMainWindow::updatePotentialImage(){
    if (checkpotentialArrayExists()){
            {
            QMutexLocker gatekeeper(&potentialLock);
            // create new empty image with appropriate dimensions
            potentialImage = QImage(potential.get_dimj(), potential.get_dimi(), QImage::Format_ARGB32);
            }

            // update sliders to match dimensions of potential, which also triggers a redraw of the image
            this->ui->slider_slicemin->setMinimum(1);
            this->ui->slider_slicemax->setMinimum(1);
            this->ui->slider_bothSlices->setMinimum(1);
            this->ui->slider_slicemin->setValue(1);
            this->ui->slider_bothSlices->setValue(1);
            this->ui->slider_slicemin->setMaximum(potential.get_dimk());
            this->ui->slider_slicemax->setMaximum(potential.get_dimk());
            this->ui->slider_bothSlices->setMaximum(potential.get_dimk());

            updatePotentialFloatImage();
        }
}

void PRISMMainWindow::updatePotentialFloatImage(){
//    if (potentialReady){
    if (checkpotentialArrayExists()){
//        QMutexLocker gatekeeper(&potentialLock);
        QMutexLocker gatekeeper(&potentialLock);

        // integrate image into the float array, then convert to uchar
        size_t min_layer = this->ui->slider_slicemin->value();
        size_t max_layer = this->ui->slider_slicemax->value();
        potentialImage_float = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{potential.get_dimj(), potential.get_dimi()}});
        for (auto k = min_layer; k <= max_layer; ++k){
            for (auto j = 0; j < potential.get_dimj(); ++j){
                for (auto i = 0; i < potential.get_dimi(); ++i){
                    potentialImage_float.at(j,i) += potential.at(k - 1,j ,i);
                }
            }
        }


        // get max/min values for contrast setting
        auto minval = std::min_element(potentialImage_float.begin(),
                                       potentialImage_float.end());
        auto maxval = std::max_element(potentialImage_float.begin(),
                                       potentialImage_float.end());
        if (ui->checkBox_sqrtIntensityPot->isChecked()){ 
            contrast_potentialMin = std::sqrt(*minval);
            contrast_potentialMax = std::sqrt(*maxval);
        } else {
            contrast_potentialMin = *minval;
            contrast_potentialMax = *maxval;
        }
        ui->lineEdit_contrastPotMin->setText(QString::number(contrast_potentialMin));
        ui->lineEdit_contrastPotMax->setText(QString::number(contrast_potentialMax));
    }
    updatePotentialDisplay();
}

void PRISMMainWindow::updatePotentialDisplay(){
//    if (potentialReady){
    if (checkpotentialArrayExists()){
            QMutexLocker gatekeeper(&potentialLock);
//            QMutexLocker gatekeeper(&dataLock);

            if (ui->checkBox_sqrtIntensityPot->isChecked()){ 
                for (auto j = 0; j < potential.get_dimj(); ++j){
                    for (auto i = 0; i < potential.get_dimi(); ++i){
//                        uchar val = getUcharFromFloat(std::sqrt(potentialImage_float.at(j,i)),
//                                                      contrast_potentialMin,
//                                                      contrast_potentialMax);
//                        potentialImage.setPixel(j, i, qRgba(val,val,val,255));
                        potentialImage.setPixel(j, i, this->colormapper.getColor(std::sqrt(potentialImage_float.at(j,i)),
                                                                                           contrast_potentialMin,
                                                                                           contrast_potentialMax));
                    }
                }
            } else {
                for (auto j = 0; j < potential.get_dimj(); ++j){
                    for (auto i = 0; i < potential.get_dimi(); ++i){
//                        uchar val = getUcharFromFloat(potentialImage_float.at(j,i),
//                                                      contrast_potentialMin,
//                                                      contrast_potentialMax);
//                        potentialImage.setPixel(j, i, qRgba(val,val,val,255));
                        potentialImage.setPixel(j, i, this->colormapper.getColor(potentialImage_float.at(j,i),
                                                                                 contrast_potentialMin,
                                                                                 contrast_potentialMax));
                    }
                }
            }


        QImage potentialImage_tmp = potentialImage.scaled(ui->lbl_image_potential->width(),
                               ui->lbl_image_potential->height(),
                               Qt::KeepAspectRatio);

        QPixmap qpix = QPixmap::fromImage( potentialImage.scaled(ui->lbl_image_potential->width(),
                                                                 ui->lbl_image_potential->height(),
                                                                 Qt::KeepAspectRatio));

        ui->lbl_image_potential_2->setPixmap(qpix);
        // draw a rectangle around the region that will be scanned
        QPainter p;
        p.begin(&qpix);
        p.setPen(QPen(Qt::yellow, 2, Qt::DotLine));
        p.drawRect(QRect(QPoint(qpix.width()  * std::max((PRISMATIC_FLOAT_PRECISION)0.0, meta->scanWindowYMin),
                                qpix.height() * std::max((PRISMATIC_FLOAT_PRECISION)0.0, meta->scanWindowXMin)),
                         QPoint(qpix.width()  * std::min((PRISMATIC_FLOAT_PRECISION)0.9999, meta->scanWindowYMax),
                                qpix.height() * std::min((PRISMATIC_FLOAT_PRECISION)0.9999, meta->scanWindowXMax))));
        p.end();
        ui->lbl_image_potential->setPixmap(qpix);


        probeImage = potentialImage;
        PRISMATIC_FLOAT_PRECISION xc, yc;
        xc = currently_calculated_X / (meta->realspacePixelSize[1]);
        yc = currently_calculated_Y / (meta->realspacePixelSize[0]);
        long xc_im, yc_im;
        xc_im = (xc / potential.get_dimi()) * probeImage.height();
        yc_im = (yc / potential.get_dimj()) * probeImage.width();

        // draw the reticle
//        const long linehalfwidth = 1;
//        const long linelength = 10;
        const long linehalfwidth = 1+probeImage.width()/200;
        const long linelength = 1+probeImage.width()/10;
        for (auto ypix = -linehalfwidth; ypix <= linehalfwidth; ++ypix){
            for (auto x = -linelength; x <= linelength; ++x){
//                std::cout << "ypix + yc_im) % probeImage.width() = " << (ypix + yc_im) % probeImage.width()<< std::endl;
             probeImage.setPixel((probeImage.width()  + (ypix + yc_im) % probeImage.width()) % probeImage.width(),
                                 (probeImage.height() + (x + xc_im) % probeImage.height()) % probeImage.height(), qRgba(255, 255, 255, 50));

            }
        }

        for (auto xpix = -linehalfwidth; xpix <= linehalfwidth; ++xpix){
            for (auto y = -linelength; y <= linelength; ++y){
             probeImage.setPixel((probeImage.width()  + (y + yc_im) % probeImage.width()) % probeImage.width(),
                                 (probeImage.height() + (xpix + xc_im) % probeImage.height()) % probeImage.height(), qRgba(255, 255, 255, 50));
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

void PRISMMainWindow::displayErrorReadingAtomsDialog(){
    QMessageBox* popup = new QMessageBox;
    popup->setWindowTitle("Prismatic: Error!");
    popup->setText(QString::fromStdString(std::string("An error occurred (see Prismatic output for more details) while attempting to read atomic coordinates from file:\n\n") +
                                          meta->filenameAtoms +
                                          std::string("\n\nEnsure that the file is accessible and is formatted correctly. Here is an example:\n\nComment line goes here\n\t5.43    5.43    5.43\n14  0.0000  0.0000  0.0000  1.0  0.076\n14  2.7150  2.7150  0.0000  1.0  0.076\n14  1.3575  4.0725  1.3575  1.0  0.076\n14  4.0725  1.3575  1.3575  1.0  0.076\n14  2.7150  0.0000  2.7150  1.0  0.076\n14  0.0000  2.7150  2.7150  1.0  0.076\n14  1.3575  1.3575  4.0725  1.0  0.076\n14  4.0725  4.0725  4.0725  1.0  0.076\n-1\n")));
    popup->show();
}

void PRISMMainWindow::displayErrorReadingParamsDialog(){
    QMessageBox* popup = new QMessageBox;
    popup->setWindowTitle("Prismatic: Error!");
    popup->setText(QString::fromStdString(std::string("An error occurred (see Prismatic output for more details) while attempting to read parameter file:\n\n") +
                                          std::string("\n\nEnsure that the file is accessible and is formatted correctly.\n")));
    popup->show();
}

void PRISMMainWindow::updateProbeK_PRISMDisplay(){
    if (checkProbesCalculated()){
        probeImage_pk = QImage(probeImage_pk_float.get_dimj(), probeImage_pk_float.get_dimi(), QImage::Format_ARGB32);
        auto contrast = std::minmax_element(probeImage_pk_float.begin(), probeImage_pk_float.end());
        double cHigh, cLow;
        cLow  = *contrast.first;
        cHigh = *contrast.second;
        if (ui->checkBox_log->isChecked()){
            cLow  = std::log(1e-5 + std::abs(cLow));
            cHigh = std::log(1e-5 + std::abs(cHigh));
            for (auto j = 0; j < probeImage_pk_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_pk_float.get_dimi(); ++i){
                    probeImage_pk.setPixel(j, i,this->colormapper.getColor(std::log(1e-5 + std::abs(probeImage_pk_float.at(j,i))),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        } else {
            for (auto j = 0; j < probeImage_pk_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_pk_float.get_dimi(); ++i){
                    probeImage_pk.setPixel(j, i,this->colormapper.getColor(probeImage_pk_float.at(j,i),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        }
        ui->lbl_image_probe_pk->setPixmap(QPixmap::fromImage(probeImage_pk.scaled(ui->lbl_image_probe_pk->width(),
                                                                                  ui->lbl_image_probe_pk->height(),
                                                                                  Qt::KeepAspectRatio)));
    }
}
void PRISMMainWindow::updateProbeR_PRISMDisplay(){
    if (checkProbesCalculated()){
        probeImage_pr = QImage(probeImage_pr_float.get_dimj(), probeImage_pr_float.get_dimi(), QImage::Format_ARGB32);
        auto contrast = std::minmax_element(probeImage_pr_float.begin(), probeImage_pr_float.end());
        double cHigh, cLow;
        cLow  = *contrast.first;
        cHigh = *contrast.second;
        if (ui->checkBox_log->isChecked()){
            cLow  = std::log(1e-5 + std::abs(cLow));
            cHigh = std::log(1e-5 + std::abs(cHigh));
            for (auto j = 0; j < probeImage_pr_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_pr_float.get_dimi(); ++i){
                    probeImage_pr.setPixel(j, i,this->colormapper.getColor(std::log(1e-5 + std::abs(probeImage_pr_float.at(j,i))),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        } else {
            for (auto j = 0; j < probeImage_pr_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_pr_float.get_dimi(); ++i){
                    probeImage_pr.setPixel(j, i,this->colormapper.getColor(probeImage_pr_float.at(j,i),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        }
        ui->lbl_image_probe_pr->setPixmap(QPixmap::fromImage(probeImage_pr.scaled(ui->lbl_image_probe_pr->width(),
                                                                                  ui->lbl_image_probe_pr->height(),
                                                                                  Qt::KeepAspectRatio)));
    }
}
void PRISMMainWindow::updateProbeK_MultisliceDisplay(){
    if (checkProbesCalculated()){
        probeImage_mk = QImage(probeImage_mk_float.get_dimj(), probeImage_mk_float.get_dimi(), QImage::Format_ARGB32);
        auto contrast = std::minmax_element(probeImage_mk_float.begin(), probeImage_mk_float.end());
        double cHigh, cLow;
        cLow  = *contrast.first;
        cHigh = *contrast.second;
        if (ui->checkBox_log->isChecked()){ 
            cLow  = std::log(1e-5 + std::abs(cLow));
            cHigh = std::log(1e-5 + std::abs(cHigh));
            for (auto j = 0; j < probeImage_mk_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_mk_float.get_dimi(); ++i){
                    probeImage_mk.setPixel(j, i,this->colormapper.getColor(std::log(1e-5 + std::abs(probeImage_mk_float.at(j,i))),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        } else {
            for (auto j = 0; j < probeImage_mk_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_mk_float.get_dimi(); ++i){
                    probeImage_mk.setPixel(j, i,this->colormapper.getColor(probeImage_mk_float.at(j,i),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        }
        ui->lbl_image_probe_mk->setPixmap(QPixmap::fromImage(probeImage_mk.scaled(ui->lbl_image_probe_mk->width(),
                                                                                  ui->lbl_image_probe_mk->height(),
                                                                                  Qt::KeepAspectRatio)));
    }
}
void PRISMMainWindow::updateProbeR_MultisliceDisplay(){
    if (checkProbesCalculated()){
        probeImage_mr = QImage(probeImage_mr_float.get_dimj(), probeImage_mr_float.get_dimi(), QImage::Format_ARGB32);
        auto contrast = std::minmax_element(probeImage_mr_float.begin(), probeImage_mr_float.end());
        double cHigh, cLow;
        cLow  = *contrast.first;
        cHigh = *contrast.second;
        if (ui->checkBox_log->isChecked()){ 
            cLow  = std::log(1e-5 + std::abs(cLow));
            cHigh = std::log(1e-5 + std::abs(cHigh));
            for (auto j = 0; j < probeImage_mr_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_mr_float.get_dimi(); ++i){
                    probeImage_mr.setPixel(j, i,this->colormapper.getColor(std::log(1e-5 + std::abs(probeImage_mr_float.at(j,i))),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        } else {
            for (auto j = 0; j < probeImage_mr_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_mr_float.get_dimi(); ++i){
                    probeImage_mr.setPixel(j, i,this->colormapper.getColor(probeImage_mr_float.at(j,i),
                                                                           cLow,
                                                                           cHigh));
                }
            }
        }
        ui->lbl_image_probe_mr->setPixmap(QPixmap::fromImage(probeImage_mr.scaled(ui->lbl_image_probe_mr->width(),
                                                                                  ui->lbl_image_probe_mr->height(),
                                                                                  Qt::KeepAspectRatio)));
    }
}

void PRISMMainWindow::updateProbe_diffRDisplay(){
    if (checkProbesCalculated()){
        probeImage_diffr = QImage(probeImage_diffr_float.get_dimj(), probeImage_diffr_float.get_dimi(), QImage::Format_ARGB32);
        auto contrast = std::minmax_element(probeImage_mr_float.begin(), probeImage_mr_float.end());
        double cHigh, cLow;
        cLow  = *contrast.first;
        cHigh = *contrast.second;
        if (ui->checkBox_log->isChecked()){
            cLow  = std::log(1e-5 + std::abs(cLow));
            cHigh = std::log(1e-5 + std::abs(cHigh));
            for (auto j = 0; j < probeImage_diffr_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_diffr_float.get_dimi(); ++i){
                    probeImage_diffr.setPixel(j, i,this->colormapper.getColor(std::log(1e-5 + std::abs(probeImage_diffr_float.at(j,i))),
                                                                              cLow,
                                                                              cHigh));
                }
            }
        } else {
            for (auto j = 0; j < probeImage_diffr_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_diffr_float.get_dimi(); ++i){
                    probeImage_diffr.setPixel(j, i,this->colormapper.getColor(probeImage_diffr_float.at(j,i),
                                                                              cLow,
                                                                              cHigh));
                }
            }
        }
        ui->lbl_image_probeDifferenceR->setPixmap(QPixmap::fromImage(probeImage_diffr.scaled(ui->lbl_image_probeDifferenceR->width(),
                                                                                             ui->lbl_image_probeDifferenceR->height(),
                                                                                             Qt::KeepAspectRatio)));
    }
}
void PRISMMainWindow::updateProbe_diffKDisplay(){
    if (checkProbesCalculated()){
        probeImage_diffk = QImage(probeImage_diffk_float.get_dimj(), probeImage_diffk_float.get_dimi(), QImage::Format_ARGB32);
        auto contrast = std::minmax_element(probeImage_mk_float.begin(), probeImage_mk_float.end());
        double cHigh, cLow;
        cLow  = *contrast.first;
        cHigh = *contrast.second;
        if (ui->checkBox_log->isChecked()){
            cLow  = std::log(1e-5 + std::abs(cLow));
            cHigh = std::log(1e-5 + std::abs(cHigh));
            for (auto j = 0; j < probeImage_diffk_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_diffk_float.get_dimi(); ++i){
                    probeImage_diffk.setPixel(j, i,this->colormapper.getColor(std::log(1e-5 + std::abs(probeImage_diffk_float.at(j,i))),
                                                                              cLow,
                                                                              cHigh));
                }
            }
        } else {
            for (auto j = 0; j < probeImage_diffk_float.get_dimj(); ++j){
                for (auto i = 0; i < probeImage_diffk_float.get_dimi(); ++i){
                    probeImage_diffk.setPixel(j, i,this->colormapper.getColor(probeImage_diffk_float.at(j,i),
                                                                              cLow,
                                                                              cHigh));
                }
            }
        }
        ui->lbl_image_probeDifferenceK->setPixmap(QPixmap::fromImage(probeImage_diffk.scaled(ui->lbl_image_probeDifferenceK->width(),
                                                                                             ui->lbl_image_probeDifferenceK->height(),
                                                                                             Qt::KeepAspectRatio)));
    }
}

void PRISMMainWindow::updateOutputImage(){
    if (checkoutputArrayExists()){
            {
            QMutexLocker gatekeeper(&outputLock);

            // create new empty image with appropriate dimensions
            outputImage = QImage(output.get_dimk(), output.get_dimj(), QImage::Format_ARGB32);
            }
            // update sliders to match dimensions of output, which also triggers a redraw of the image
            this->ui->slider_angmin->setMinimum(0);
            this->ui->slider_angmax->setMinimum(0);
            this->ui->slider_bothSlices->setMinimum(0);
            this->ui->slider_angmin->setMaximum(detectorAngles.size() - 1);
            this->ui->slider_angmax->setMaximum(detectorAngles.size() - 1);
            this->ui->slider_bothSlices->setMaximum(detectorAngles.size() - 1);
            this->ui->slider_angmax->setValue(std::min(this->ui->slider_angmax->value(), this->ui->slider_angmax->maximum()));
            this->ui->lineEdit_angmin->setText(QString::number(std::max((PRISMATIC_FLOAT_PRECISION)0.0, detectorAngles[this->ui->slider_angmin->value()] - (detectorAngles[1] - detectorAngles[0])/2)));
            this->ui->lineEdit_angmax->setText(QString::number(detectorAngles[std::min((int)detectorAngles.size(), this->ui->slider_angmax->value())] + (detectorAngles[1] - detectorAngles[0])/2));
        }
    updateOutputFloatImage();
}

void PRISMMainWindow::updateOutputFloatImage(){
      if (checkoutputArrayExists()){
        QMutexLocker gatekeeper(&outputLock);

        // integrate image into the float array, then convert to uchar
        size_t min_layer = this->ui->slider_angmin->value();
        size_t max_layer = this->ui->slider_angmax->value();
        outputImage_float = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{output.get_dimk(), output.get_dimj()}});
        for (auto j = 0; j < output.get_dimk(); ++j){
            for (auto i = 0; i < output.get_dimj(); ++i){
                 for (auto k = min_layer; k <= max_layer; ++k){
                    outputImage_float.at(j,i) += output.at(0, j, i, k);
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
        ui->lineEdit_contrast_outputMin->setText(QString::number(contrast_outputMin));
        ui->lineEdit_contrast_outputMax->setText(QString::number(contrast_outputMax));
    }
    updateOutputDisplay();
}

void PRISMMainWindow::updateOutputDisplay(){
//    if (outputReady){
    if (checkoutputArrayExists()){
        QMutexLocker gatekeeper(&outputLock);
            for (auto j = 0; j < output.get_dimk(); ++j){
                for (auto i = 0; i < output.get_dimj(); ++i){
//                    uchar val = getUcharFromFloat(outputImage_float.at(j,i),
//                                                  contrast_outputMin,
//                                                  contrast_outputMax);
//                    outputImage.setPixel(j, i, qRgba(val,val,val,255));
                    outputImage.setPixel(j, i,this->colormapper.getColor(outputImage_float.at(j,i),
                                                                         contrast_outputMin,
                                                                         contrast_outputMax));

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

void PRISMMainWindow::updateOutputImage_HRTEM(){
    if (checkoutputArrayExists_HRTEM()){
            {
            QMutexLocker gatekeeper(&outputLock);

            // create new empty image with appropriate dimensions
            outputImage_HRTEM = QImage(smatrix.get_dimj(), smatrix.get_dimi(), QImage::Format_ARGB32);
            }
            // update sliders to match dimensions of output, which also triggers a redraw of the image
            this->ui->slider_angmin_2->setMinimum(0);
            this->ui->slider_angmin_2->setMaximum(smatrix.get_dimk()-1);
            this->ui->lineEdit_angmin_2->setText(QString::number((PRISMATIC_FLOAT_PRECISION) smatrix.get_dimk()));
        }
    updateOutputFloatImage_HRTEM();
}

void PRISMMainWindow::updateOutputFloatImage_HRTEM(){
      if (checkoutputArrayExists_HRTEM()){
        QMutexLocker gatekeeper(&outputLock);

        // integrate image into the float array, then convert to uchar
        size_t beam = 0; // this->ui->slider_angmin_2->value();
        outputImage_HRTEM_float = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{smatrix.get_dimj(), smatrix.get_dimi()}});
        for (auto j = 0; j < smatrix.get_dimj(); ++j){
            for (auto i = 0; i < smatrix.get_dimi(); ++i){
                outputImage_HRTEM_float.at(j,i) += pow(std::abs(smatrix.at(beam, j, i)), 2.0);
            }
        }

        // get max/min values for contrast setting
        auto minval = std::min_element(outputImage_HRTEM_float.begin(),
                                       outputImage_HRTEM_float.end());
        auto maxval = std::max_element(outputImage_HRTEM_float.begin(),
                                       outputImage_HRTEM_float.end());
        contrast_outputMin_HRTEM = *minval;
        contrast_outputMax_HRTEM = *maxval;
        ui->lineEdit_contrast_outputMin_2->setText(QString::number(contrast_outputMin));
        ui->lineEdit_contrast_outputMax_2->setText(QString::number(contrast_outputMax));
    }
    updateOutputDisplay_HRTEM();
}

void PRISMMainWindow::updateOutputDisplay_HRTEM(){
    if (checkoutputArrayExists_HRTEM()){
        QMutexLocker gatekeeper(&outputLock);
            for (auto j = 0; j < smatrix.get_dimj(); ++j){
                for (auto i = 0; i < smatrix.get_dimi(); ++i){
                    outputImage_HRTEM.setPixel(j, i, this->colormapper.getColor(outputImage_HRTEM_float.at(j,i),
                                                                         contrast_outputMin_HRTEM,
                                                                         contrast_outputMax_HRTEM));
                }
            }

        QImage outputImage_tmp = outputImage_HRTEM.scaled(ui->lbl_image_output_2->width(),
                                                    ui->lbl_image_output_2->height(),
                                                    Qt::KeepAspectRatio);

        ui->lbl_image_output_2->setPixmap(QPixmap::fromImage( outputImage_HRTEM.scaled(ui->lbl_image_output_2->width(),
                                                                               ui->lbl_image_output_2->height(),
                                                                               Qt::KeepAspectRatio)));
    }
}

void PRISMMainWindow::updateProbeImages(){
    updateProbeR_PRISMDisplay();
    updateProbeK_PRISMDisplay();
    updateProbeR_MultisliceDisplay();
    updateProbeK_MultisliceDisplay();
    updateProbe_diffRDisplay();
    updateProbe_diffKDisplay();
}

void PRISMMainWindow::updateAllImages(){
    updateOutputImage();
    updatePotentialImage();
    updateProbeImages();
}

void PRISMMainWindow::updateSliders_fromLineEdits(){
    this->ui->slider_slicemin->setValue(std::min(this->ui->lineEdit_slicemin->text().toInt(),
                                                 this->ui->slider_slicemax->value()));
    this->ui->slider_bothSlices->setValue(this->ui->slider_slicemin->value());
    this->ui->slider_slicemax->setValue(std::max(this->ui->lineEdit_slicemax->text().toInt(),
                                                 this->ui->slider_slicemin->value()));
}
void PRISMMainWindow::updateSliders_fromLineEdits_ang(){
//    if (outputReady){
    if (checkoutputArrayExists()){
        bool flagMin = false;
        bool flagMax = false;
        PRISMATIC_FLOAT_PRECISION minval = ( (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_angmin->text().toDouble(&flagMin)) /
        (detectorAngles[1]-detectorAngles[0]);
        PRISMATIC_FLOAT_PRECISION maxval = ( (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_angmax->text().toDouble(&flagMax)) /
        (detectorAngles[1]-detectorAngles[0]);
        if (flagMin & flagMax){
            this->ui->slider_angmin->setValue(std::min( (int)std::round(minval),
                                                         this->ui->slider_angmax->value()));
            this->ui->slider_angmax->setValue(std::max( (int)std::round(maxval) - 1,
                                                         this->ui->slider_angmin->value()));
            this->ui->slider_bothDetectors->setValue(this->ui->slider_angmin->value());
            updateSlider_lineEdits_max_ang(this->ui->slider_angmax->value());
            updateSlider_lineEdits_min_ang(this->ui->slider_angmin->value());
        }
    }
}

void PRISMMainWindow::updateSlider_lineEdits_min(int val){
    if (val <= this->ui->slider_slicemax->value()){
        this->ui->lineEdit_slicemin->setText(QString::number(val));
    } else {
        this->ui->slider_slicemin->setValue(this->ui->slider_slicemax->value());
    }
    this->ui->slider_bothSlices->setValue(this->ui->slider_slicemin->value());

}

void PRISMMainWindow::updateSlider_lineEdits_max(int val){
    if (val >= this->ui->slider_slicemin->value()){
        this->ui->lineEdit_slicemax->setText(QString::number(val));
    } else {
        this->ui->slider_slicemax->setValue(this->ui->slider_slicemin->value());
    }
}

void PRISMMainWindow::updateSlider_lineEdits_max_ang(int val){

    if (checkoutputArrayExists()){
        if (val >= this->ui->slider_angmin->value()){
            double scaled_val = detectorAngles[0] + val * (detectorAngles[1] - detectorAngles[0]) + (detectorAngles[1] - detectorAngles[0])/2;
            this->ui->lineEdit_angmax->setText(QString::number(scaled_val));
        } else {
            this->ui->slider_angmax->setValue(this->ui->slider_angmin->value());
        }        
    }
}

void PRISMMainWindow::updateSlider_lineEdits_min_ang(int val){
    if (checkoutputArrayExists()){
        if (val <= this->ui->slider_angmax->value()){
            double scaled_val = detectorAngles[0] + val * (detectorAngles[1] - detectorAngles[0]) - (detectorAngles[1] - detectorAngles[0])/2;
            this->ui->lineEdit_angmin->setText(QString::number(scaled_val));
        } else {
            this->ui->slider_angmin->setValue(this->ui->slider_angmax->value());
        }
        this->ui->slider_bothDetectors->setValue(this->ui->slider_angmin->value());
    }
}


void PRISMMainWindow::updateContrastPotMin(){
    bool flag = false;
    contrast_potentialMin = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_contrastPotMin->text().toDouble(&flag);
    if (flag)updatePotentialDisplay();
}
void PRISMMainWindow::updateContrastPotMax(){
    bool flag = false;
    contrast_potentialMax = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_contrastPotMax->text().toDouble(&flag);
    if (flag)updatePotentialDisplay();
}

void PRISMMainWindow::updateContrastAngMin(){
    bool flag = false;
    contrast_outputMin = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_contrast_outputMin->text().toDouble(&flag);
    if (flag)updateOutputDisplay();
}
void PRISMMainWindow::updateContrastAngMax(){
    bool flag = false;
    contrast_outputMax = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_contrast_outputMax->text().toDouble(&flag);
    if (flag)updateOutputDisplay();
}

void PRISMMainWindow::updateAlphaMax(){
    using namespace Prismatic;
    PRISMATIC_FLOAT_PRECISION f_x = 4 * meta->interpolationFactorX;
    PRISMATIC_FLOAT_PRECISION f_y = 4 * meta->interpolationFactorY;
    Array1D<size_t> imageSize({{(size_t)(meta->cellDim[1] * meta->tileY), (size_t)(meta->cellDim[2] * meta->tileX)}}, {{2}});
    imageSize[0] = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)4.0,  (PRISMATIC_FLOAT_PRECISION)(f_y * round(((PRISMATIC_FLOAT_PRECISION)imageSize[0]) / meta->realspacePixelSize[0] / f_y)));
    imageSize[1] = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)4.0,  (PRISMATIC_FLOAT_PRECISION)(f_x * round(((PRISMATIC_FLOAT_PRECISION)imageSize[1]) / meta->realspacePixelSize[1] / f_x)));

    long long ncx = (size_t) floor((PRISMATIC_FLOAT_PRECISION) imageSize[1] / 2);
    PRISMATIC_FLOAT_PRECISION dpx = 1.0 / ((PRISMATIC_FLOAT_PRECISION)imageSize[1] * meta->realspacePixelSize[1]);
    long long ncy = (size_t) floor((PRISMATIC_FLOAT_PRECISION) imageSize[0] / 2);
    PRISMATIC_FLOAT_PRECISION dpy = 1.0 / ((PRISMATIC_FLOAT_PRECISION)imageSize[0] * meta->realspacePixelSize[0]);
    PRISMATIC_FLOAT_PRECISION qMax = std::min(dpx*(ncx), dpy*(ncy)) / 2;

    PRISMATIC_FLOAT_PRECISION alphaMax = qMax * calculateLambda(*meta);
    ui->lbl_alphaMax->setText(QString::fromUtf8("\u03B1 max = ") + QString::number(alphaMax * 1000) + QString(" mrad"));
}

void PRISMMainWindow::checkInput_lineEdit_scanWindowXMin(){
    bool flagMin = false;
    bool flagMax = false;
    PRISMATIC_FLOAT_PRECISION minVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowXMin->text().toDouble(&flagMin);
    PRISMATIC_FLOAT_PRECISION maxVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowXMax->text().toDouble(&flagMax);
    if (flagMin & flagMax){
        if (minVal > maxVal){
            meta->scanWindowXMin = meta->scanWindowXMax;
            ui->lineEdit_scanWindowXMin->setText(QString::number(meta->scanWindowXMin));
        }
    } else {
        meta->scanWindowXMin = 0.0;
        ui->lineEdit_scanWindowXMin->setText(QString::number(meta->scanWindowXMin));
    }
    if (!minWindowYSet){
        this->ui->lineEdit_scanWindowYMin->setText(QString::number(meta->scanWindowXMin));
        this->meta->scanWindowYMin = meta->scanWindowXMin;
    }
}

void PRISMMainWindow::checkInput_lineEdit_scanWindowXMax(){
    bool flagMin = false;
    bool flagMax = false;
    PRISMATIC_FLOAT_PRECISION minVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowXMin->text().toDouble(&flagMin);
    PRISMATIC_FLOAT_PRECISION maxVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowXMax->text().toDouble(&flagMax);
    if (flagMin & flagMax){
        if (maxVal < minVal){
            meta->scanWindowXMax = meta->scanWindowXMin;
            ui->lineEdit_scanWindowXMax->setText(QString::number(meta->scanWindowXMax));
        }
    } else {
        meta->scanWindowXMax = 1.0;
        ui->lineEdit_scanWindowXMax->setText(QString::number(meta->scanWindowXMax));
    }
    if (!maxWindowYSet){
        this->ui->lineEdit_scanWindowYMax->setText(QString::number(meta->scanWindowXMax));
        this->meta->scanWindowYMax = meta->scanWindowXMax;
    }
}

void PRISMMainWindow::checkInput_lineEdit_scanWindowYMin(){
    bool flagMin = false;
    bool flagMax = false;
    PRISMATIC_FLOAT_PRECISION minVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowYMin->text().toDouble(&flagMin);
    PRISMATIC_FLOAT_PRECISION maxVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowYMax->text().toDouble(&flagMax);
    if (flagMin & flagMax){
        if (minVal > maxVal){
            meta->scanWindowYMin = meta->scanWindowYMax;
            ui->lineEdit_scanWindowYMin->setText(QString::number(meta->scanWindowYMin));
        }
    } else {
        meta->scanWindowYMin = 0.0;
        ui->lineEdit_scanWindowYMin->setText(QString::number(meta->scanWindowYMin));
    }
}

void PRISMMainWindow::checkInput_lineEdit_scanWindowYMax(){
    bool flagMin = false;
    bool flagMax = false;
    PRISMATIC_FLOAT_PRECISION minVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowYMin->text().toDouble(&flagMin);
    PRISMATIC_FLOAT_PRECISION maxVal = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_scanWindowYMax->text().toDouble(&flagMax);
    if (flagMin & flagMax){
        if (maxVal < minVal){
            meta->scanWindowYMax = meta->scanWindowYMin;
            ui->lineEdit_scanWindowYMax->setText(QString::number(meta->scanWindowYMax));
        }
    } else {
        meta->scanWindowYMax = 1.0;
        ui->lineEdit_scanWindowYMax->setText(QString::number(meta->scanWindowYMax));
    }
}
void PRISMMainWindow::checkInput_lineEdit_cellDimX(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_cellDimX->text().toDouble(&flag);
    if (!flag | (val < 0)){
        meta->cellDim[2] = 1;
        ui->lineEdit_cellDimX->setText(QString::number(meta->cellDim[2]));
    }
}

void PRISMMainWindow::checkInput_lineEdit_cellDimY(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_cellDimY->text().toDouble(&flag);
    if (!flag | (val < 0)){
        meta->cellDim[1] = 1;
        ui->lineEdit_cellDimY->setText(QString::number(meta->cellDim[1]));
    }
}

void PRISMMainWindow::checkInput_lineEdit_cellDimZ(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_cellDimZ->text().toDouble(&flag);
    if (!flag | (val < 0)){
        meta->cellDim[0] = 1;
        ui->lineEdit_cellDimZ->setText(QString::number(meta->cellDim[0]));
    }
}

void PRISMMainWindow::checkInput_lineEdit_tileX(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_tileX->text().toInt(&flag);
    if (!flag | (val < 1)){
        meta->tileX = 1;
    }
    ui->lineEdit_tileX->setText(QString::number(meta->tileX));
}

void PRISMMainWindow::checkInput_lineEdit_tileY(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_tileY->text().toInt(&flag);
    if (!flag | (val < 1)){
        meta->tileY = 1;
    }
    ui->lineEdit_tileY->setText(QString::number(meta->tileY));
}

void PRISMMainWindow::checkInput_lineEdit_tileZ(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_tileZ->text().toInt(&flag);
    if (!flag | (val < 1)){
        meta->tileZ = 1;
    }
    ui->lineEdit_tileZ->setText(QString::number(meta->tileZ));
}


void PRISMMainWindow::checkInput_lineEdit_interpFactor_x(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_interpFactor_x->text().toInt(&flag);
    if (!flag | (val < 1)){
        meta->interpolationFactorX = 1;
    }
    ui->lineEdit_interpFactor_x->setText(QString::number(meta->interpolationFactorX));
    if (!interpYSet){
        this->ui->lineEdit_interpFactor_y->setText(QString::number(meta->interpolationFactorX));
        this->meta->interpolationFactorY = meta->interpolationFactorX;
    }
}

void PRISMMainWindow::checkInput_lineEdit_interpFactor_y(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_interpFactor_y->text().toInt(&flag);
    if (!flag | (val < 1)){
        meta->interpolationFactorY = 1;
    }
    ui->lineEdit_interpFactor_y->setText(QString::number(meta->interpolationFactorY));
}
void PRISMMainWindow::checkInput_lineEdit_pixelSizeX(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_pixelSizeX->text().toDouble(&flag);
    if (!flag | (val < 0)){
        meta->realspacePixelSize[1] = 0.1;
    }
    ui->lineEdit_pixelSizeX->setText(QString::number(meta->realspacePixelSize[1]));
    if (!pixelSizeYSet){
        this->ui->lineEdit_pixelSizeY->setText(QString::number(meta->realspacePixelSize[1]));
        this->meta->realspacePixelSize[0] = meta->realspacePixelSize[1];
    }
}


void PRISMMainWindow::checkInput_lineEdit_pixelSizeY(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)ui->lineEdit_pixelSizeY->text().toDouble(&flag);
    if (!flag | (val < 0)){
        meta->realspacePixelSize[0] = 0.1;
    }
    ui->lineEdit_pixelSizeY->setText(QString::number(meta->realspacePixelSize[0]));
}
void PRISMMainWindow::saveCurrentOutputImage(){
    if (checkoutputArrayExists()){
        QMutexLocker gatekeeper(&outputLock);
        uint8_t * tmp_buffer = (uint8_t*) malloc(outputImage_float.size()*8);
        PRISMATIC_FLOAT_PRECISION max_val = *std::max_element(outputImage_float.begin(),outputImage_float.end());
        PRISMATIC_FLOAT_PRECISION min_val = *std::min_element(outputImage_float.begin(),outputImage_float.end());
        for( auto i = 0; i < outputImage_float.size(); i++) tmp_buffer[i] = std::floor(255*(outputImage_float[i]-min_val)/(max_val-min_val));
        uchar * buffer = (uchar*) &tmp_buffer[0];
        QImage image(buffer,outputImage_float.get_dimi(),outputImage_float.get_dimj(),outputImage_float.get_dimi()*1,QImage::Format_Grayscale8);
        image.save(ui->lineEdit_saveOutputImage->text().toStdString().c_str(),"PNG",100);
        free(tmp_buffer);  
        //outputImage_float.toMRC_f(ui->lineEdit_saveOutputImage->text().toStdString().c_str());
    }
}

void PRISMMainWindow::setStreamingMode(int val){
    enum{Auto=0, SingleXfer=1, Stream=2} setting;
    switch (val){
        case Auto:
            meta->transferMode = Prismatic::StreamingMode::Auto;
            std::cout << "Setting streaming mode: Auto" << std::endl;
            break;
        case SingleXfer:
            meta->transferMode = Prismatic::StreamingMode::SingleXfer;
            std::cout << "Setting streaming mode: Single Transfer" << std::endl;
            break;
        case Stream:
            meta->transferMode = Prismatic::StreamingMode::Stream;
            std::cout << "Setting streaming mode: Streaming" << std::endl;
            break;
    }
}

void PRISMMainWindow::newRandomSeed(){
    PRISMATIC_FLOAT_PRECISION val = rand() % 100000;
    ui->lineEdit_randomSeed->setText(QString::number(val));
    meta->randomSeed = val;
}

void PRISMMainWindow::toggle2DOutput(){
    meta->save2DOutput = ui->checkBox_2D->isChecked();
    if(meta->save2DOutput){
        ui->lineEdit_2D_inner->setEnabled(true);
        ui->lineEdit_2D_outer->setEnabled(true);
    }else{
        ui->lineEdit_2D_inner->setDisabled(true);
        ui->lineEdit_2D_outer->setDisabled(true);
    }
}

void PRISMMainWindow::toggle3DOutput(){
    meta->save3DOutput = ui->checkBox_3D->isChecked();
}

void PRISMMainWindow::toggle4DOutput(){
    meta->save4DOutput = ui->checkBox_4D->isChecked();
}

void PRISMMainWindow::toggle4Dcrop(){
    meta->crop4DOutput = ui->checkBox_4D->isChecked();
    if(meta->crop4DOutput){
        ui->lineEdit_4Damax->setEnabled(true);
    }else{
        ui->lineEdit_4Damax->setDisabled(true);
    }
}

void PRISMMainWindow::set4Damax_fromLineEdit(){
    bool flag = false;
    PRISMATIC_FLOAT_PRECISION val = (PRISMATIC_FLOAT_PRECISION)this->ui->lineEdit_4Damax->text().toDouble(&flag);
    if (flag){
        this->meta->crop4Damax =  val / 1000;
        std::cout << "Setting 4D crop angle to " << val << " mrad" << std::endl;
    }
    resetCalculation();
}

void PRISMMainWindow::toggleDPC_CoM(){
    meta->saveDPC_CoM = ui->checkBox_DPC_CoM->isChecked();
}

void PRISMMainWindow::togglePotentialSlices(){
    meta->savePotentialSlices = ui->checkBox_PS->isChecked();
}

void PRISMMainWindow::toggleSMatrixoutput(){
    meta->saveSMatrix = ui->checkBox_saveSMatrix->isChecked();
}

void PRISMMainWindow::toggleComplexoutput(){
    meta->saveComplexOutputWave = ui->checkBox_saveComplex->isChecked();
}

void PRISMMainWindow::toggleProbeOutput(){
    meta->saveProbe = ui->checkBox_saveProbe->isChecked();
}

void PRISMMainWindow::toggleThermalEffects(){
    meta->includeThermalEffects = ui->checkBox_thermalEffects->isChecked();
    resetCalculation();
}

void PRISMMainWindow::togglematrixRefocus(){
    meta->matrixRefocus = ui->checkBox_matrixRefocus->isChecked();
    resetCalculation();
}

void PRISMMainWindow::togglePotential3D(){
    meta->potential3D = ui->checkBox_potential3D->isChecked();
    resetCalculation();
}

/*void PRISMMainWindow::toggleOccupancy(){
    meta->includeOccupancy = ui->checkBox_occupancy->isChecked();
    resetCalculation();
}*/

void PRISMMainWindow::toggleNyquist(){
    meta->nyquistSampling = ui->checkBox_NQS->isChecked();
    if(meta->nyquistSampling){
        ui->lineEdit_probeStepX->setDisabled(true);
        ui->lineEdit_probeStepY->setDisabled(true);
    }else{
        ui->lineEdit_probeStepX->setEnabled(true);
        ui->lineEdit_probeStepY->setEnabled(true);
    }
    resetCalculation();
}

void PRISMMainWindow::changeColormap(QString color){
    this->colormapper.setColormap(Prismatic::ColorTable[color]);
//    this->colormapper.setColormap(Prismatic::ColorTable[QString("Jet")]);

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

//void PRISMMainWindow::toggleSaveProjectedPotential(){
//    this->saveProjectedPotential = ui->checkBox_saveProjectedPotential->isChecked() ? true:false;
//}

bool PRISMMainWindow::potentialIsReady(){
    QMutexLocker gatekeeper(&dataLock);
    return potentialReady;
}

bool PRISMMainWindow::overwriteFile(){
    QMutexLocker gatekeeper(&dataLock);
    return overwriteCheck;
}

bool PRISMMainWindow::SMatrixIsReady(){
//    return false;
    QMutexLocker gatekeeper(&dataLock);
    return ScompactReady;
}
bool PRISMMainWindow::OutputIsReady(){
    QMutexLocker gatekeeper(&dataLock);
    return outputReady;
}

bool PRISMMainWindow::checkoutputArrayExists(){
    QMutexLocker gatekeeper(&outputLock);
    return outputArrayExists;
}

bool PRISMMainWindow::checkoutputArrayExists_HRTEM(){
    QMutexLocker gatekeeper(&outputLock);
    return outputArrayExists_HRTEM;
}

bool PRISMMainWindow::checkpotentialArrayExists(){
    QMutexLocker gatekeeper(&potentialLock);
    return potentialArrayExists;
}

bool PRISMMainWindow::checkProbesCalculated(){
    QMutexLocker gatekeeper(&probeLock);
    return probesCalculated;
}

void PRISMMainWindow::probeK_PRISMReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr){
    {
        QMutexLocker gatekeeper(&probeLock);
        probeImage_pk_float = arr;
        probesCalculated = true;
    }
    updateProbeK_PRISMDisplay();
}

void PRISMMainWindow::probeR_PRISMReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr){
    {
        QMutexLocker gatekeeper(&probeLock);
        probeImage_pr_float = arr;
        probesCalculated = true;
    }
    updateProbeR_PRISMDisplay();
}

void PRISMMainWindow::probeK_MultisliceReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr){
    {
        QMutexLocker gatekeeper(&probeLock);
        probeImage_mk_float = arr;
        probesCalculated = true;
    }
    updateProbeK_MultisliceDisplay();
}

void PRISMMainWindow::probeR_MultisliceReceived(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr){
    {
        QMutexLocker gatekeeper(&probeLock);
        probeImage_mr_float = arr;
        probesCalculated = true;
    }
    updateProbeR_MultisliceDisplay();
}

void PRISMMainWindow::calculateProbe_diffR(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr){
    {
        QMutexLocker gatekeeper(&probeLock);
        probeImage_diffr_float = arr;
        probesCalculated = true;
    }
    updateProbe_diffRDisplay();
}

void PRISMMainWindow::calculateProbe_diffK(Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> arr){
    {
        QMutexLocker gatekeeper(&probeLock);
        probeImage_diffk_float = arr;
        probesCalculated = true;
    }
    updateProbe_diffKDisplay();
}

void PRISMMainWindow::flipOverwrite(){
    {
        QMutexLocker gatekeeper(&dataLock);
        if(overwriteCheck){
            overwriteCheck = false;
        }else{
            overwriteCheck = true;
        }
    }
}

void PRISMMainWindow::potentialReceived(Prismatic::Array3D<PRISMATIC_FLOAT_PRECISION> _potential){
    {
        QMutexLocker gatekeeper(&potentialLock);
        potential = _potential;
        potentialArrayExists = true;
    }
    {
        QMutexLocker gatekeeper(&dataLock);
        potentialReady = true;
    }
}

void PRISMMainWindow::outputReceived(Prismatic::Array4D<PRISMATIC_FLOAT_PRECISION> _output){
    {
        QMutexLocker gatekeeper(&outputLock);
        output = _output;
        outputArrayExists = true;
    }
}

void PRISMMainWindow::outputReceived_HRTEM(Prismatic::Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> _output){
    {
        QMutexLocker gatekeeper(&outputLock);
        smatrix = _output;
        outputArrayExists_HRTEM = true;
    }
}

void PRISMMainWindow::moveBothPotentialSliders(int val){
    int difference = ui->slider_slicemax->value() - ui->slider_slicemin->value();
    if (val + difference <= ui->slider_slicemax->maximum()){
        ui->slider_slicemax->setValue(val + difference);
        ui->slider_slicemin->setValue(val);
    }
    ui->slider_bothSlices->setValue(ui->slider_slicemin->value());
}

void PRISMMainWindow::moveBothDetectorSliders(int val){
    int difference = ui->slider_angmax->value() - ui->slider_angmin->value();
    if (val + difference <= ui->slider_angmax->maximum()){
        ui->slider_angmax->setValue(val + difference);
        ui->slider_angmin->setValue(val);
    }
    ui->slider_bothDetectors->setValue(ui->slider_angmin->value());
}


void PRISMMainWindow::enableOutputWidgets(){
    ui->slider_angmax->setEnabled(true);
    ui->slider_angmin->setEnabled(true);
    ui->slider_bothDetectors->setEnabled(true);
    ui->lineEdit_contrast_outputMax->setEnabled(true);
    ui->lineEdit_contrast_outputMin->setEnabled(true);
    ui->lineEdit_angmin->setEnabled(true);
    ui->lineEdit_angmax->setEnabled(true);
}
void PRISMMainWindow::resetCalculation(){
//    QMutexLocker gatekeeper2(&calculationLock);
    QMutexLocker gatekeeper(&dataLock);
    std::cout << "Resetting Calculation" << std::endl;
    outputReady     = false;
    ScompactReady   = false;
    potentialReady  = false;
    probeSetupReady = false;
    // delete this->meta;
    // this->meta = new Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION>;
    // updateDisplay();
}

void PRISMMainWindow::resetPotential(){
    QMutexLocker gatekeeper(&dataLock);
    potentialReady  = false;
}


void PRISMMainWindow::resetLinks(){
    // reset linked parameters
    minWindowYSet   = false;
    maxWindowYSet   = false;
    interpYSet      = false;
    pixelSizeYSet   = false;
    probeStepYSet   = false;
    probeTiltYSet   = false;
}

void PRISMMainWindow::redrawImages(){

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

    ui->lbl_image_output_2->setPixmap(QPixmap::fromImage(outputImage_HRTEM.scaled(ui->lbl_image_output_2->width(),
                                                                          ui->lbl_image_output_2->height(),
                                                                          Qt::KeepAspectRatio)));

    updatePotentialDisplay();
}

void PRISMMainWindow::setscan_WindowYMin_edited(){
    minWindowYSet = true;
}

void PRISMMainWindow::setscan_WindowYMax_edited(){
    maxWindowYSet = true;
}

void PRISMMainWindow::setinterpYSet_edited(){
    interpYSet = true;
}

void PRISMMainWindow::setpixelSizeYSet_edited(){
    pixelSizeYSet = true;
}

void PRISMMainWindow::setprobeStepYSet_edited(){
    probeStepYSet = true;
}

void PRISMMainWindow::setprobeTiltYSet_edited(){
    probeTiltYSet = true;
}

void PRISMMainWindow::userHasSetCellDims(){
    meta->userSpecifiedCelldims = true;
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

unsigned char getUcharFromFloat(PRISMATIC_FLOAT_PRECISION val,
                                PRISMATIC_FLOAT_PRECISION contrast_low,
                                PRISMATIC_FLOAT_PRECISION contrast_high){

    if (val < contrast_low)  return 0;
    if (val > contrast_high) return 255;
    return (unsigned char)( (val - contrast_low) / (contrast_high - contrast_low) * 255);

}

PRISMATIC_FLOAT_PRECISION calculateLambda(Prismatic::Metadata<PRISMATIC_FLOAT_PRECISION> meta){
	constexpr double m = 9.109383e-31;
	constexpr double e = 1.602177e-19;
	constexpr double c = 299792458;
	constexpr double h = 6.62607e-34;
	return (PRISMATIC_FLOAT_PRECISION)(h / sqrt(2 * m * e * meta.E0) / sqrt(1 + e * meta.E0 / 2 / m / c / c) * 1e10);
}


//Collapses and Opens Sample Box window
void PRISMMainWindow::collapseSample(){

    //if widget is open
    if(!sampleClosed){
        this->ui->box_samplesettings->setMaximumHeight(boxClosed);
        this->ui->box_samplesettings->setMinimumHeight(boxClosed);
        this->ui->btn_closebox->setText("+");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_samplesettings, "minimumHeight" );
        //pAni->setStartValue(boxOpen);
        //pAni->setEndValue(boxClosed);
        //pAni->setDuration(animSpeed);
        //pAni->start();

        sampleClosed = true;


    }else{

        this->ui->box_samplesettings->setMaximumHeight(boxOpen);
        this->ui->box_samplesettings->setMinimumHeight(boxOpen);
        this->ui->btn_closebox->setText("-");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_samplesettings, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(boxOpen);
        //pAni->setDuration(animSpeed);

        //QPropertyAnimation* pAni2 = new QPropertyAnimation(this->ui->scrollArea_4, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(scrollOpen);
        //pAni->setDuration(300);
        this->ui->scrollArea_4->setMaximumHeight(scrollOpen);
        this->ui->scrollArea_4->setMinimumHeight(scrollOpen);

        //pAni->start();
        //pAni2->start();
        sampleClosed = false;
    }
}


void PRISMMainWindow::collapseStem(){

    //if widget is open
    if(!stemClosed){
        this->ui->box_stemsettings->setMaximumHeight(boxClosed);
        this->ui->box_stemsettings->setMinimumHeight(boxClosed);
        this->ui->btn_closebox_3->setText("+");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_stemsettings, "minimumHeight" );
        //pAni->setStartValue(boxOpen);
        //pAni->setEndValue(boxClosed);
        //pAni->setDuration(animSpeed);
        //pAni->start();
        
        stemClosed = true;


    }else{

        this->ui->box_stemsettings->setMaximumHeight(boxOpen);
        this->ui->box_stemsettings->setMinimumHeight(boxOpen);
        this->ui->btn_closebox_3->setText("-");


        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_stemsettings, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(boxOpen);
        //pAni->setDuration(animSpeed);

        this->ui->scrollArea_3->setMinimumHeight(scrollOpen);
        this->ui->scrollArea_3->setMaximumHeight(scrollOpen);

        //pAni->start();
        stemClosed = false;
    }
}

void PRISMMainWindow::collapseHrtem(){

    //if widget is open
    if(!hrtemClosed){
        this->ui->box_hrtemsettings->setMaximumHeight(boxClosed);
        this->ui->box_hrtemsettings->setMinimumHeight(boxClosed);
        this->ui->btn_closebox_4->setText("+");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_hrtemsettings, "minimumHeight" );
        //pAni->setStartValue(boxOpen);
        //pAni->setEndValue(boxClosed);
        //pAni->setDuration(animSpeed);
        //pAni->start();
        
        hrtemClosed = true;


    }else{
        this->ui->box_hrtemsettings->setMaximumHeight(boxOpen);
        this->ui->box_hrtemsettings->setMinimumHeight(boxOpen);
        this->ui->btn_closebox_4->setText("-");


        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_hrtemsettings, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(boxOpen);
        //pAni->setDuration(animSpeed);

        this->ui->scrollArea_7->setMinimumHeight(scrollOpen);
        this->ui->scrollArea_7->setMaximumHeight(scrollOpen);

        //pAni->start();
        hrtemClosed = false;
    }
}


void PRISMMainWindow::collapseOutput(){

    //if widget is open
    if(!outputClosed){
        this->ui->box_outputsettings->setMaximumHeight(boxClosed);
        this->ui->box_outputsettings->setMinimumHeight(boxClosed);
        this->ui->btn_closebox_5->setText("+");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_outputsettings, "minimumHeight" );
        //pAni->setStartValue(120);
        //pAni->setEndValue(boxClosed);
        //pAni->setDuration(animSpeed);
        //pAni->start();

        outputClosed = true;


    }else{

        this->ui->box_outputsettings->setMaximumHeight(120);
        this->ui->box_outputsettings->setMinimumHeight(120);
        this->ui->btn_closebox_5->setText("-");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_outputsettings, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(120);
        //pAni->setDuration(animSpeed);

        this->ui->scrollArea_6->setMinimumHeight(110);
        this->ui->scrollArea_6->setMaximumHeight(110);

        //pAni->start();

        outputClosed = false;
    }
}

void PRISMMainWindow::collapseSimulation(){

    //if widget is open
    if(!simulationClosed){


        this->ui->box_simulationsettings->setMaximumHeight(boxClosed);
        this->ui->box_simulationsettings->setMinimumHeight(boxClosed);
        this->ui->btn_closebox_2->setText("+");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_simulationsettings, "minimumHeight" );
        //pAni->setStartValue(boxOpen);
        //pAni->setEndValue(boxClosed);
        //pAni->setDuration(animSpeed);
        //pAni->start();

        simulationClosed = true;


    }else{


        this->ui->box_simulationsettings->setMaximumHeight(boxOpen+10);
        this->ui->box_simulationsettings->setMinimumHeight(boxOpen+10);
        this->ui->btn_closebox_2->setText("-");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_simulationsettings, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(boxOpen);
        //pAni->setDuration(animSpeed);

        this->ui->scrollArea_5->setMaximumHeight(scrollOpen+5);
        this->ui->scrollArea_5->setMinimumHeight(scrollOpen+5);

        //pAni->start();

        simulationClosed = false;
    }
}

void PRISMMainWindow::collapseComputational(){

    //if widget is open
    if(!computationalClosed){

        this->ui->box_computationalsettings->setMaximumHeight(boxClosed);
        this->ui->box_computationalsettings->setMinimumHeight(boxClosed);

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_computationalsettings, "minimumHeight" );
        //pAni->setStartValue(190);
        //pAni->setEndValue(boxClosed);
        //pAni->setDuration(animSpeed);
        //pAni->start();

        computationalClosed = true;
        this->ui->btn_closebox_7->setText("+");


    }else{

        this->ui->box_computationalsettings->setMaximumHeight(200);
        this->ui->box_computationalsettings->setMinimumHeight(200);
        this->ui->btn_closebox_7->setText("-");

        //QPropertyAnimation* pAni = new QPropertyAnimation(this->ui->box_computationalsettings, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(190);
        //pAni->setDuration(animSpeed);


        //QPropertyAnimation* pAni2 = new QPropertyAnimation(this->ui->scrollArea_9, "minimumHeight" );
        //pAni->setStartValue(boxClosed);
        //pAni->setEndValue(170);
        //pAni->setDuration(300);


        this->ui->scrollArea_9->setMaximumHeight(170);
        this->ui->scrollArea_9->setMinimumHeight(170);


        //pAni->start();
        //pAni2->start();

        computationalClosed = false;
    }
}


void PRISMMainWindow::lightField(){

// set stylesheet
    QFile file(":/light.qss");
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream stream(&file);
    qApp->setStyleSheet(stream.readAll());

}

void PRISMMainWindow::darkField(){

// set stylesheet
    QFile file(":/dark.qss");
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream stream(&file);
    qApp->setStyleSheet(stream.readAll());

}
