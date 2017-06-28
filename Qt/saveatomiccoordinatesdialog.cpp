#include "saveatomiccoordinatesdialog.h"
#include "ui_saveatomiccoordinatesdialog.h"

SaveAtomicCoordinatesDialog::SaveAtomicCoordinatesDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SaveAtomicCoordinatesDialog)
{
    ui->setupUi(this);
}

void SaveAtomicCoordinatesDialog::SaveAtomCoords(){
    emit signalSaveAtomCoords(ui->lineEdit_Filename->text(), ui->lineEdit_Comment->text());
}

SaveAtomicCoordinatesDialog::~SaveAtomicCoordinatesDialog()
{
    delete ui;
}
