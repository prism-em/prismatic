#ifndef PRISMMAINWINDOW_H
#define PRISMMAINWINDOW_H

#include <QMainWindow>
#include <iostream>
#include "meta.h"

namespace Ui {
	class PRISMMainWindow;
}

class PRISMMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit PRISMMainWindow(QWidget *parent = 0);
    ~PRISMMainWindow();
public slots:
	void setInterpolationFactor();
	void setFilenameAtoms_fromLineEdit();
	void setFilenameAtoms_fromDialog();
	void setFilenameOutput_fromLineEdit();
	void setFilenameOutput_fromDialog();
	void setFilenameAtoms(const std::string& filename);
	void setFilenameOutput(const std::string& filename);
	void launch();
private:
    Ui::PRISMMainWindow *ui;
    PRISM::Metadata<double> *meta;

};

#endif // PRISMMAINWINDOW_H
