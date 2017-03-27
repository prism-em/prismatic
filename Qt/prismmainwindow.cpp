#include <QFileDialog>
#include "prismmainwindow.h"
#include "ui_prismmainwindow.h"
#include <fstream>
#include <iostream>
//#include "PRISM_entry.h"
#include "configure.h"
#include "prism_qthreads.h"

PRISMMainWindow::PRISMMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PRISMMainWindow),
    potentialReady(false),
    potentialImage(QImage())
{
	// build Qt generated interface
    ui->setupUi(this);

	// set window title
	this->setWindowTitle("PRISM");


    ui->box_sampleSettings->setStyleSheet("QGroupBox { \
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
                                         left: 150px;\
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

    QPixmap potentialImage("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/Qt/potential.png");
    QPixmap probeImage("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/Qt/probe.png");
    QPixmap outputImage("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/Qt/output.png");
    potentialScene = new QGraphicsScene(this);
    QGraphicsScene* probeScene = new QGraphicsScene(this);
    QGraphicsScene* outputScene = new QGraphicsScene(this);
    potentialImage = potentialImage.scaled(ui->graphicsView_potential->width()*0.9, ui->graphicsView_potential->height()*0.9, Qt::KeepAspectRatio);
    //probeImage     = probeImage.scaled(ui->graphicsView_probe->width()*0.9, ui->graphicsView_probe->height()*0.9, Qt::KeepAspectRatio);
   probeImage     = probeImage.scaled(ui->graphicsView_potential->width()*0.9, ui->graphicsView_potential->height()*0.9, Qt::KeepAspectRatio);

    outputImage    = outputImage.scaled(ui->graphicsView_potential->width()*0.9, ui->graphicsView_potential->height()*0.9, Qt::KeepAspectRatio);



   // potentialScene->addPixmap(potentialImage);
    probeScene->addPixmap(probeImage);
    outputScene->addPixmap(outputImage);


    //potentialScene->setSceneRect(potentialImage.rect());
    probeScene->setSceneRect(probeImage.rect());
    outputScene->setSceneRect(outputImage.rect());

//    potentialScene->setSceneRect(0,0,250,250);
   // probeScene->setSceneRect(probeImage.rect());
  //  outputScene->setSceneRect(outputImage.rect());
    ui->graphicsView_potential->setScene(potentialScene);
   // ui->graphicsView_potential->fitInView(potentialImage.rect());
    ui->graphicsView_probe->setScene(probeScene);
   // ui->graphicsView_probe->setScene(potentialScene);
    ui->graphicsView_output->setScene(outputScene);

    int s = 512;
    QImage image(s, s, QImage::Format_ARGB32);
    for (int i = 0; i < s; ++i){
        for (int j = 0; j < s; ++j){
            image.setPixel(i, j,  qRgba( i, j, 0, 255));
        }
    }
    potentialScene->addPixmap(QPixmap::fromImage(image));
   // potentialScene->setSceneRect(image.rect());

	// set initially displayed values based on the default parameters
	this->meta = new PRISM::Metadata<PRISM_FLOAT_PRECISION>;
	{
		std::stringstream ss;
        ss << this->meta->interpolationFactor;
        this->ui->lineedit_interpFactor_x->setText(QString::fromStdString(ss.str()));
		ss.str("");
		ss << this->meta->potBound;
		this->ui->lineEdit_potbound->setText(QString::fromStdString(ss.str()));
		ss.str("");
		ss << this->meta->alphaBeamMax;
		this->ui->lineEdit_alphaBeamMax->setText(QString::fromStdString(ss.str()));
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
		ss << (this->meta->E0 *1e-3);
		this->ui->lineEdit_E0->setText(QString::fromStdString(ss.str()));
		ss.str("");

		this->ui->lineedit_outputfile->setText(QString::fromStdString(ss.str()));
		this->ui->spinBox_numGPUs->setValue(this->meta->NUM_GPUS);
		this->ui->spinBox_numThreads->setValue(this->meta->NUM_THREADS);
		this->ui->spinBox_numFP->setValue(this->meta->numFP);

	}

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
#endif //PRISM_ENABLE_GPU

	this->ui->lineedit_outputfile->setText(QString::fromStdString(this->meta->filename_output));
    connect(this->ui->lineedit_interpFactor_x,SIGNAL(editingFinished()),this,SLOT(setInterpolationFactor()));
	connect(this->ui->lineedit_outputfile,SIGNAL(editingFinished()),this,SLOT(setFilenameOutput_fromLineEdit()));
	connect(this->ui->btn_atomsfile_browse, SIGNAL(pressed()), this, SLOT(setFilenameAtoms_fromDialog()));
    //connect(this->ui->btn_outputfile_browse, SIGNAL(pressed()), this, SLOT(setFilenameOutput_fromDialog()));
	connect(this->ui->btn_go,SIGNAL(pressed()),this,SLOT(launch()));
    connect(this->ui->spinBox_numGPUs, SIGNAL(valueChanged(int)), this, SLOT(setNumGPUs(const int&)));
    connect(this->ui->spinBox_numThreads, SIGNAL(valueChanged(int)), this, SLOT(setNumThreads(const int&)));
    connect(this->ui->spinBox_numFP, SIGNAL(valueChanged(int)), this, SLOT(setNumFP(const int&)));
    connect(this->ui->lineEdit_alphaBeamMax, SIGNAL(editingFinished()), this, SLOT(setAlphaBeamMax_fromLineEdit()));
    connect(this->ui->lineedit_pixelSize, SIGNAL(editingFinished()), this, SLOT(setPixelSize_fromLineEdit()));
    connect(this->ui->lineEdit_potbound, SIGNAL(editingFinished()), this, SLOT(setPotBound_fromLineEdit()));
    connect(this->ui->lineEdit_sliceThickness, SIGNAL(editingFinished()), this, SLOT(setSliceThickness_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimX, SIGNAL(editingFinished()), this, SLOT(setCellDimX_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimY, SIGNAL(editingFinished()), this, SLOT(setCellDimY_fromLineEdit()));
    connect(this->ui->lineEdit_cellDimZ, SIGNAL(editingFinished()), this, SLOT(setCellDimZ_fromLineEdit()));
    connect(this->ui->lineEdit_E0, SIGNAL(editingFinished()), this, SLOT(setE0_fromLineEdit()));
	connect(this->ui->radBtn_PRISM, SIGNAL(clicked(bool)), this, SLOT(setAlgo_PRISM()));
	connect(this->ui->radBtn_Multislice, SIGNAL(clicked(bool)), this, SLOT(setAlgo_Multislice()));
    connect(this->ui->btn_calcPotential, SIGNAL(clicked(bool)), this, SLOT(calculatePotential()));
    connect(this->ui->lineEdit_slicemin, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits()));
    connect(this->ui->lineEdit_slicemax, SIGNAL(editingFinished()), this, SLOT(updateSliders_fromLineEdits()));
    connect(this->ui->slider_slicemin, SIGNAL(valueChanged(int)), this, SLOT(updateSlider_lineEdits_min(int)));
    connect(this->ui->slider_slicemax, SIGNAL(valueChanged(int)), this, SLOT(updateSlider_lineEdits_max(int)));


//    connect(this->ui->btn_calcSmatrix, SIGNAL(clicked(bool)), this, SLOT(testImage()));

}

//    void PRISMMainWindow::testImage(){
//        static int offset = 0;
//        std::cout << "offset = " << offset << std::endl;

//        int s = 1024;
//        QImage image(s, s, QImage::Format_ARGB32);
//        for (int i = 0; i < s; ++i){
//            for (int j = 0; j < s; ++j){
//                image.setPixel(i, j,  qRgba(i + offset, j + offset, offset, 255));
//            }
//        }
//        offset++;
//        this->potentialScene->addPixmap(QPixmap::fromImage(image));
//    }


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
}

void PRISMMainWindow::setInterpolationFactor(){
	bool flag;
    const size_t& new_f = this->ui->lineedit_interpFactor_x->text().toUInt(&flag);
	if (flag){
		std::cout << "Setting interpolation factor to " << new_f<< std::endl;
		this->meta->interpolationFactor = new_f;
	} else{
        std::cout << "Invalid interpolation factor input: " <<  this->ui->lineedit_interpFactor_x->text().toStdString() << std::endl;
	}
}

void PRISMMainWindow::setFilenameAtoms_fromDialog(){
	QString filename;
	filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename);
	this->setFilenameAtoms(filename.toStdString());
}

void PRISMMainWindow::setFilenameOutput_fromDialog(){
	QString filename;
	filename = QFileDialog::getOpenFileName(this, tr("AnyFile"), filename);
	this->ui->lineedit_outputfile->setText(filename);
	this->setFilenameOutput(filename.toStdString());
}


void PRISMMainWindow::setFilenameOutput_fromLineEdit(){
	const std::string& filename = this->ui->lineedit_outputfile->text().toStdString();
	this->setFilenameOutput(filename);
}

void PRISMMainWindow::setFilenameAtoms(const std::string& filename){
	std::cout << "Setting atoms filename to " << filename << std::endl;
	this->meta->filename_atoms = filename;
}

void PRISMMainWindow::setFilenameOutput(const std::string& filename){
	std::cout << "Setting output filename to " << filename << std::endl;
	this->meta->filename_output = filename;
}

void PRISMMainWindow::launch(){
	std::cout << "Launching PRISM calculation with the following paramters:\n";
	std::cout << "Atoms filename = " << this->meta->filename_atoms << '\n';
	std::cout << "Output filename = " << this->meta->filename_output << '\n';
	std::cout << "Interpolation factor = " << this->meta->interpolationFactor << '\n';
	std::cout << "pixelSize[0] = " << this->meta->realspace_pixelSize<< '\n';

    PRISM::configure(*this->meta);
	int returnCode =  PRISM::execute_plan(*this->meta);
	if (returnCode == 0){
		std::cout << "Calculation complete" << std::endl;
	} else {
		std::cout << "Calculation returned error code " << returnCode << std::endl;
	}

}

void PRISMMainWindow::setNumGPUs(const int& num){
    if (num > 0){
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

}

void PRISMMainWindow::setPixelSize_fromLineEdit(){
    PRISM_FLOAT_PRECISION val =(PRISM_FLOAT_PRECISION)this->ui->lineedit_pixelSize->text().toDouble();
    if (val > 0){
        this->meta->realspace_pixelSize = val;
        this->meta->realspace_pixelSize = val;
        std::cout << "Setting X/Y pixel size to " << val << std::endl;
    }
}

void PRISMMainWindow::setPotBound_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_potbound->text().toDouble();
    if (val > 0){
        this->meta->potBound = val;
        std::cout << "Setting potential bound to " << val << std::endl;
    }
}

void PRISMMainWindow::setAlphaBeamMax_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_alphaBeamMax->text().toDouble();
    if (val > 0){
        this->meta->alphaBeamMax = val;
        std::cout << "Setting alphaBeamMax to " << val << std::endl;
    }
}

void PRISMMainWindow::setSliceThickness_fromLineEdit(){
    PRISM_FLOAT_PRECISION val = (PRISM_FLOAT_PRECISION)this->ui->lineEdit_sliceThickness->text().toDouble();
    if (val > 0){
        this->meta->sliceThickness = val;
        std::cout << "Setting sliceThickness to " << val << std::endl;
    }
}

void PRISMMainWindow::setCellDimX_fromLineEdit(){
    int val = this->ui->lineEdit_cellDimX->text().toInt();
    if (val > 0){
        this->meta->cellDim[2] = (size_t)val;
        std::cout << "Setting X cell dimension to " << val << std::endl;
    }
}

void PRISMMainWindow::setCellDimY_fromLineEdit(){
    int val = this->ui->lineEdit_cellDimY->text().toInt();
    if (val > 0){
        this->meta->cellDim[1] = (size_t)val;
        std::cout << "Setting Y cell dimension to " << val << std::endl;
    }
}

void PRISMMainWindow::setCellDimZ_fromLineEdit(){
    int val = this->ui->lineEdit_cellDimZ->text().toInt();
    if (val > 0){
        this->meta->cellDim[0] = (size_t)val;
        std::cout << "Setting Z cell dimension to " << val << std::endl;
    }
}

void PRISMMainWindow::setE0_fromLineEdit(){
    int val = this->ui->lineEdit_E0->text().toInt();
    if (val > 0){
        this->meta->E0 = val * 1e3;
        std::cout << "Setting E0 to " << val << std::endl;
    }
}

void PRISMMainWindow::calculatePotential(){
PotentialThread *worker = new PotentialThread(this);
std::cout <<"starting working" << std::endl;
worker->meta.toString();
worker->start();
std::cout <<"worker started" << std::endl;
connect(worker, SIGNAL(finished()), this, SLOT(updatePotentialImage()));
connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));

}

void PRISMMainWindow::updatePotentialImage(){
    std::cout << "updatePotentialImage called" << std::endl;
    if (potentialReady){
        potentialImage = QImage(potential.get_dimj(), potential.get_dimi(), QImage::Format_Grayscale8);
        this->ui->slider_slicemin->setMinimum(0);
        this->ui->slider_slicemax->setMinimum(0);
        this->ui->slider_slicemin->setMaximum(potential.get_dimk());
        this->ui->slider_slicemax->setMaximum(potential.get_dimk());
        this->ui->slider_slicemax->setValue(potential.get_dimk());
    }
}

void PRISMMainWindow::updateSliders_fromLineEdits(){
    this->ui->slider_slicemin->setValue(std::min(this->ui->lineEdit_slicemin->text().toInt(),
                                                 this->ui->slider_slicemax->value()));
    this->ui->slider_slicemax->setValue(std::max(this->ui->lineEdit_slicemax->text().toInt(),
                                                 this->ui->slider_slicemin->value()));
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

PRISMMainWindow::~PRISMMainWindow()
{
    delete ui;
	delete meta;
}
