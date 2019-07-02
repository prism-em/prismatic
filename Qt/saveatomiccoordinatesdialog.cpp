#include "saveatomiccoordinatesdialog.h"
#include "ui_saveatomiccoordinatesdialog.h"

SaveAtomicCoordinatesDialog::SaveAtomicCoordinatesDialog(QWidget *parent) : QDialog(parent),
                                                                            ui(new Ui::SaveAtomicCoordinatesDialog)
{
    ui->setupUi(this);
}

void SaveAtomicCoordinatesDialog::SaveAtomCoords()
{
    emit signalSaveAtomCoords(ui->lineEdit_Filename->text(), ui->lineEdit_Comment->text());
}

void SaveAtomicCoordinatesDialog::setFilenameText(QString text)
{
    ui->lineEdit_Filename->setText(text);
}

void SaveAtomicCoordinatesDialog::setCommentText(QString text)
{
    ui->lineEdit_Comment->setText(text);
}

SaveAtomicCoordinatesDialog::~SaveAtomicCoordinatesDialog()
{
    delete ui;
}
