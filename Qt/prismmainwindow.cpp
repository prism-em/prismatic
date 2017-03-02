#include <QFileDialog>
#include "prismmainwindow.h"
#include "ui_prismmainwindow.h"
#include <fstream>
#include <iostream>
#include "PRISM_entry.h"

PRISMMainWindow::PRISMMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PRISMMainWindow)
{
    ui->setupUi(this);
	this->meta = new PRISM::Metadata<double>;

    this->setWindowTitle("PRISM");
	connect(this->ui->lineedit_f,SIGNAL(editingFinished()),this,SLOT(setInterpolationFactor()));
	connect(this->ui->lineedit_atomsfile,SIGNAL(editingFinished()),this,SLOT(setFilenameAtoms_fromLineEdit()));
	connect(this->ui->lineedit_outputfile,SIGNAL(editingFinished()),this,SLOT(setFilenameOutput_fromLineEdit()));
	connect(this->ui->btn_atomsfile_browse, SIGNAL(pressed()), this, SLOT(setFilenameAtoms_fromDialog()));
	connect(this->ui->btn_outputfile_browse, SIGNAL(pressed()), this, SLOT(setFilenameOutput_fromDialog()));
	connect(this->ui->btn_go,SIGNAL(pressed()),this,SLOT(launch()));
}

void PRISMMainWindow::setInterpolationFactor(){
	bool flag;
	const size_t& new_f = this->ui->lineedit_f->text().toUInt(&flag);
	if (flag){
		std::cout << "Setting interpolation factor to " << new_f<< std::endl;
		this->meta->interpolationFactor = new_f;
	} else{
		std::cout << "Invalid interpolation factor input: " <<  this->ui->lineedit_f->text().toStdString() << std::endl;
	}
}

void PRISMMainWindow::setFilenameAtoms_fromDialog(){
	QString filename;
	filename = QFileDialog::getOpenFileName(this, tr("ExistingFile"), filename);
	this->ui->lineedit_atomsfile->text(filename);
	this->setFilenameAtoms(filename.toStdString());
}

void PRISMMainWindow::setFilenameOutput_fromDialog(){
	QString filename;
	filename = QFileDialog::getOpenFileName(this, tr("AnyFile"), filename);
	this->ui->lineedit_outputfile->text(filename);
	this->setFilenameOutput(filename.toStdString());
}

void PRISMMainWindow::setFilenameAtoms_fromLineEdit(){
	const std::string& filename = this->ui->lineedit_atomsfile->text().toStdString();
	this->setFilenameAtoms(filename);
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

	PRISM::PRISM_entry((*this->meta));
}



PRISMMainWindow::~PRISMMainWindow()
{
    delete ui;
	delete meta;
}
