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
	void launch();
    void setNumGPUs(const int& numGPUs);
    void setNumThreads(const int& numThreads);
    void setNumFP(const int& numFP);
    void setPixelSize_fromLineEdit();
    void setPotBound_fromLineEdit();
    void setAlphaBeamMax_fromLineEdit();
    void setSliceThickness_fromLineEdit();
    void setCellDimX_fromLineEdit();
    void setCellDimY_fromLineEdit();
    void setCellDimZ_fromLineEdit();
    void setE0_fromLineEdit();

private:
    Ui::PRISMMainWindow *ui;
    PRISM::Metadata<double> *meta;
    void setFilenameAtoms(const std::string& filename);
    void setFilenameOutput(const std::string& filename);
    void setRealspacePixelSize(const double& pixel_size);
    void setPotBound(const double& potBound);
    void setNumFP(const size_t& numFP);
    void setE0(const double& E0);
    void setAlphaBeamMax(const double& alphaBeamMax);
    void setSliceThickness(const double& thickness);
    void setCellDimX(const int& dimX);
    void setCellDimY(const int& dimY);
    void setCellDimZ(const int& dimZ);
};

#endif // PRISMMAINWINDOW_H
