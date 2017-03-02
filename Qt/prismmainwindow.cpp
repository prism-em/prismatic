#include "prismmainwindow.h"
#include "ui_prismmainwindow.h"

PRISMMainWindow::PRISMMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PRISMMainWindow)
{
    ui->setupUi(this);
    this->setWindowTitle("PRISM");
    connect(this->ui->btn_go,SIGNAL(pressed()),this,SLOT(printMessage()));
}

PRISMMainWindow::~PRISMMainWindow()
{
    delete ui;
}
