#ifndef PRISM_PROGRESSBAR_H
#define PRISM_PROGRESSBAR_H

#include <QDialog>
#include "prismmainwindow.h"
#include "meta.h"
#include "configure.h"
namespace Ui {
class prism_progressbar;
}

class prism_progressbar : public QDialog
{
    Q_OBJECT

public:
    explicit prism_progressbar(PRISMMainWindow *parent);
    ~prism_progressbar();

private:
    Ui::prism_progressbar *ui;
    PRISMMainWindow *parent;
    prism_progressbar *progressbar;
    PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
};

#endif // PRISM_PROGRESSBAR_H
