#include "prism_progressbar.h"
#include "ui_prism_progressbar.h"

prism_progressbar::prism_progressbar(PRISMMainWindow *_parent) :
    parent(_parent),
    ui(new Ui::prism_progressbar)
{
    ui->setupUi(this);
}

prism_progressbar::~prism_progressbar()
{
    delete ui;
}
