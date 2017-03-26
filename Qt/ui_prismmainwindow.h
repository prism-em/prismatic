/********************************************************************************
** Form generated from reading UI file 'prismmainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PRISMMAINWINDOW_H
#define UI_PRISMMAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PRISMMainWindow
{
public:
    QWidget *centralWidget;
    QGroupBox *box_simulationSettings;
    QGroupBox *box_calculationSettings;
    QGroupBox *box_sampleSettings;
    QHBoxLayout *horizontalLayout_3;
    QVBoxLayout *verticalLayout_6;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *btn_loadCoordinates;
    QPushButton *btn_saveCoordinates;
    QHBoxLayout *horizontalLayout;
    QVBoxLayout *verticalLayout_5;
    QLabel *label_8;
    QLabel *label_cellDim;
    QLabel *label_6;
    QVBoxLayout *verticalLayout;
    QLabel *label_2;
    QLineEdit *lineEdit_cellDimX;
    QLineEdit *lineEdit_tileX;
    QVBoxLayout *verticalLayout_2;
    QLabel *label_3;
    QLineEdit *lineEdit_cellDimY;
    QLineEdit *lineEdit_tileY;
    QVBoxLayout *verticalLayout_3;
    QLabel *label_4;
    QLineEdit *lineEdit_cellDimZ;
    QLineEdit *lineEdit_tileZ;
    QVBoxLayout *verticalLayout_4;
    QLabel *label_9;
    QLabel *label_5;
    QLabel *label_7;
    QLabel *lbl_algo;
    QRadioButton *radBtn_PRISM;
    QRadioButton *radBtn_Multislice;
    QPushButton *btn_go;
    QLineEdit *lineedit_outputfile;
    QPushButton *btn_outputfile_browse;
    QLabel *lbl_outputfile;
    QPushButton *btn_atomsfile_browse;
    QLineEdit *lineedit_atomsfile;
    QLabel *lbl_atomsfile;
    QLabel *label;
    QLineEdit *lineedit_f;
    QLabel *lbl_pixelSize;
    QLineEdit *lineedit_pixelSize;
    QLineEdit *lineEdit_alphaBeamMax;
    QLabel *label_alphaBeamMax;
    QSpinBox *spinBox_numFP;
    QLabel *lbl_numfp;
    QSpinBox *spinBox_numGPUs;
    QLabel *lbl_numgpus;
    QLineEdit *lineEdit_potbound;
    QLabel *lbl_potBound;
    QLabel *label_sliceThickness;
    QLineEdit *lineEdit_sliceThickness;
    QSpinBox *spinBox_numThreads;
    QLabel *lbl_numthreads;
    QLabel *lbl_E0;
    QLineEdit *lineEdit_E0;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PRISMMainWindow)
    {
        if (PRISMMainWindow->objectName().isEmpty())
            PRISMMainWindow->setObjectName(QStringLiteral("PRISMMainWindow"));
        PRISMMainWindow->resize(901, 584);
        centralWidget = new QWidget(PRISMMainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        box_simulationSettings = new QGroupBox(centralWidget);
        box_simulationSettings->setObjectName(QStringLiteral("box_simulationSettings"));
        box_simulationSettings->setGeometry(QRect(560, 170, 221, 151));
        box_calculationSettings = new QGroupBox(centralWidget);
        box_calculationSettings->setObjectName(QStringLiteral("box_calculationSettings"));
        box_calculationSettings->setGeometry(QRect(550, 310, 261, 171));
        box_sampleSettings = new QGroupBox(centralWidget);
        box_sampleSettings->setObjectName(QStringLiteral("box_sampleSettings"));
        box_sampleSettings->setGeometry(QRect(550, 60, 251, 111));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(box_sampleSettings->sizePolicy().hasHeightForWidth());
        box_sampleSettings->setSizePolicy(sizePolicy);
        horizontalLayout_3 = new QHBoxLayout(box_sampleSettings);
        horizontalLayout_3->setSpacing(0);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(2, 2, 2, 2);
        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setSpacing(0);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        verticalLayout_6->setContentsMargins(0, -1, -1, -1);
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(0);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        btn_loadCoordinates = new QPushButton(box_sampleSettings);
        btn_loadCoordinates->setObjectName(QStringLiteral("btn_loadCoordinates"));

        horizontalLayout_2->addWidget(btn_loadCoordinates);

        btn_saveCoordinates = new QPushButton(box_sampleSettings);
        btn_saveCoordinates->setObjectName(QStringLiteral("btn_saveCoordinates"));

        horizontalLayout_2->addWidget(btn_saveCoordinates);


        verticalLayout_6->addLayout(horizontalLayout_2);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(0);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(0);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        label_8 = new QLabel(box_sampleSettings);
        label_8->setObjectName(QStringLiteral("label_8"));

        verticalLayout_5->addWidget(label_8);

        label_cellDim = new QLabel(box_sampleSettings);
        label_cellDim->setObjectName(QStringLiteral("label_cellDim"));
        label_cellDim->setAlignment(Qt::AlignCenter);

        verticalLayout_5->addWidget(label_cellDim);

        label_6 = new QLabel(box_sampleSettings);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setAlignment(Qt::AlignCenter);

        verticalLayout_5->addWidget(label_6);


        horizontalLayout->addLayout(verticalLayout_5);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(0);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        label_2 = new QLabel(box_sampleSettings);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setMaximumSize(QSize(16777215, 30));
        label_2->setAlignment(Qt::AlignCenter);

        verticalLayout->addWidget(label_2);

        lineEdit_cellDimX = new QLineEdit(box_sampleSettings);
        lineEdit_cellDimX->setObjectName(QStringLiteral("lineEdit_cellDimX"));
        QSizePolicy sizePolicy1(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(lineEdit_cellDimX->sizePolicy().hasHeightForWidth());
        lineEdit_cellDimX->setSizePolicy(sizePolicy1);
        lineEdit_cellDimX->setMaximumSize(QSize(40, 16777215));
        lineEdit_cellDimX->setLayoutDirection(Qt::LeftToRight);
        lineEdit_cellDimX->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout->addWidget(lineEdit_cellDimX);

        lineEdit_tileX = new QLineEdit(box_sampleSettings);
        lineEdit_tileX->setObjectName(QStringLiteral("lineEdit_tileX"));
        lineEdit_tileX->setMaximumSize(QSize(40, 16777215));
        lineEdit_tileX->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout->addWidget(lineEdit_tileX);


        horizontalLayout->addLayout(verticalLayout);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(0);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        label_3 = new QLabel(box_sampleSettings);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setMaximumSize(QSize(16777215, 30));
        label_3->setAlignment(Qt::AlignCenter);

        verticalLayout_2->addWidget(label_3);

        lineEdit_cellDimY = new QLineEdit(box_sampleSettings);
        lineEdit_cellDimY->setObjectName(QStringLiteral("lineEdit_cellDimY"));
        sizePolicy1.setHeightForWidth(lineEdit_cellDimY->sizePolicy().hasHeightForWidth());
        lineEdit_cellDimY->setSizePolicy(sizePolicy1);
        lineEdit_cellDimY->setMaximumSize(QSize(40, 16777215));
        lineEdit_cellDimY->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_2->addWidget(lineEdit_cellDimY);

        lineEdit_tileY = new QLineEdit(box_sampleSettings);
        lineEdit_tileY->setObjectName(QStringLiteral("lineEdit_tileY"));
        lineEdit_tileY->setMaximumSize(QSize(40, 16777215));
        lineEdit_tileY->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_2->addWidget(lineEdit_tileY);


        horizontalLayout->addLayout(verticalLayout_2);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(0);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        label_4 = new QLabel(box_sampleSettings);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setMaximumSize(QSize(16777215, 30));
        label_4->setAlignment(Qt::AlignCenter);

        verticalLayout_3->addWidget(label_4);

        lineEdit_cellDimZ = new QLineEdit(box_sampleSettings);
        lineEdit_cellDimZ->setObjectName(QStringLiteral("lineEdit_cellDimZ"));
        sizePolicy1.setHeightForWidth(lineEdit_cellDimZ->sizePolicy().hasHeightForWidth());
        lineEdit_cellDimZ->setSizePolicy(sizePolicy1);
        lineEdit_cellDimZ->setMaximumSize(QSize(40, 16777215));
        lineEdit_cellDimZ->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_3->addWidget(lineEdit_cellDimZ);

        lineEdit_tileZ = new QLineEdit(box_sampleSettings);
        lineEdit_tileZ->setObjectName(QStringLiteral("lineEdit_tileZ"));
        lineEdit_tileZ->setMaximumSize(QSize(40, 16777215));
        lineEdit_tileZ->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_3->addWidget(lineEdit_tileZ);


        horizontalLayout->addLayout(verticalLayout_3);

        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(0);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        label_9 = new QLabel(box_sampleSettings);
        label_9->setObjectName(QStringLiteral("label_9"));

        verticalLayout_4->addWidget(label_9);

        label_5 = new QLabel(box_sampleSettings);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setAlignment(Qt::AlignCenter);

        verticalLayout_4->addWidget(label_5);

        label_7 = new QLabel(box_sampleSettings);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setAlignment(Qt::AlignCenter);

        verticalLayout_4->addWidget(label_7);


        horizontalLayout->addLayout(verticalLayout_4);


        verticalLayout_6->addLayout(horizontalLayout);


        horizontalLayout_3->addLayout(verticalLayout_6);

        label_2->raise();
        label_3->raise();
        label_4->raise();
        label_5->raise();
        label_6->raise();
        lineEdit_tileX->raise();
        lineEdit_tileY->raise();
        lineEdit_tileZ->raise();
        label_7->raise();
        lineEdit_tileX->raise();
        label_cellDim->raise();
        label_8->raise();
        label_9->raise();
        label_9->raise();
        btn_loadCoordinates->raise();
        btn_saveCoordinates->raise();
        lbl_algo = new QLabel(centralWidget);
        lbl_algo->setObjectName(QStringLiteral("lbl_algo"));
        lbl_algo->setGeometry(QRect(13, 325, 59, 16));
        radBtn_PRISM = new QRadioButton(centralWidget);
        radBtn_PRISM->setObjectName(QStringLiteral("radBtn_PRISM"));
        radBtn_PRISM->setGeometry(QRect(78, 324, 65, 20));
        radBtn_Multislice = new QRadioButton(centralWidget);
        radBtn_Multislice->setObjectName(QStringLiteral("radBtn_Multislice"));
        radBtn_Multislice->setGeometry(QRect(153, 324, 83, 20));
        btn_go = new QPushButton(centralWidget);
        btn_go->setObjectName(QStringLiteral("btn_go"));
        btn_go->setGeometry(QRect(360, 10, 185, 145));
        btn_go->setMinimumSize(QSize(0, 145));
        QFont font;
        font.setPointSize(30);
        font.setBold(true);
        font.setWeight(75);
        btn_go->setFont(font);
        lineedit_outputfile = new QLineEdit(centralWidget);
        lineedit_outputfile->setObjectName(QStringLiteral("lineedit_outputfile"));
        lineedit_outputfile->setGeometry(QRect(89, 60, 125, 21));
        btn_outputfile_browse = new QPushButton(centralWidget);
        btn_outputfile_browse->setObjectName(QStringLiteral("btn_outputfile_browse"));
        btn_outputfile_browse->setGeometry(QRect(216, 56, 87, 32));
        lbl_outputfile = new QLabel(centralWidget);
        lbl_outputfile->setObjectName(QStringLiteral("lbl_outputfile"));
        lbl_outputfile->setGeometry(QRect(14, 60, 67, 16));
        btn_atomsfile_browse = new QPushButton(centralWidget);
        btn_atomsfile_browse->setObjectName(QStringLiteral("btn_atomsfile_browse"));
        btn_atomsfile_browse->setGeometry(QRect(216, 13, 87, 32));
        lineedit_atomsfile = new QLineEdit(centralWidget);
        lineedit_atomsfile->setObjectName(QStringLiteral("lineedit_atomsfile"));
        lineedit_atomsfile->setGeometry(QRect(89, 17, 125, 21));
        lbl_atomsfile = new QLabel(centralWidget);
        lbl_atomsfile->setObjectName(QStringLiteral("lbl_atomsfile"));
        lbl_atomsfile->setGeometry(QRect(14, 17, 67, 16));
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(14, 100, 134, 16));
        lineedit_f = new QLineEdit(centralWidget);
        lineedit_f->setObjectName(QStringLiteral("lineedit_f"));
        lineedit_f->setGeometry(QRect(156, 100, 50, 21));
        lineedit_f->setMaximumSize(QSize(50, 16777215));
        lbl_pixelSize = new QLabel(centralWidget);
        lbl_pixelSize->setObjectName(QStringLiteral("lbl_pixelSize"));
        lbl_pixelSize->setGeometry(QRect(17, 139, 58, 16));
        lineedit_pixelSize = new QLineEdit(centralWidget);
        lineedit_pixelSize->setObjectName(QStringLiteral("lineedit_pixelSize"));
        lineedit_pixelSize->setGeometry(QRect(83, 139, 125, 21));
        sizePolicy1.setHeightForWidth(lineedit_pixelSize->sizePolicy().hasHeightForWidth());
        lineedit_pixelSize->setSizePolicy(sizePolicy1);
        lineedit_pixelSize->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lineEdit_alphaBeamMax = new QLineEdit(centralWidget);
        lineEdit_alphaBeamMax->setObjectName(QStringLiteral("lineEdit_alphaBeamMax"));
        lineEdit_alphaBeamMax->setGeometry(QRect(132, 208, 125, 21));
        sizePolicy1.setHeightForWidth(lineEdit_alphaBeamMax->sizePolicy().hasHeightForWidth());
        lineEdit_alphaBeamMax->setSizePolicy(sizePolicy1);
        lineEdit_alphaBeamMax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_alphaBeamMax = new QLabel(centralWidget);
        label_alphaBeamMax->setObjectName(QStringLiteral("label_alphaBeamMax"));
        label_alphaBeamMax->setGeometry(QRect(17, 208, 107, 16));
        spinBox_numFP = new QSpinBox(centralWidget);
        spinBox_numFP->setObjectName(QStringLiteral("spinBox_numFP"));
        spinBox_numFP->setGeometry(QRect(393, 139, 48, 24));
        spinBox_numFP->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lbl_numfp = new QLabel(centralWidget);
        lbl_numfp->setObjectName(QStringLiteral("lbl_numfp"));
        lbl_numfp->setGeometry(QRect(269, 139, 116, 24));
        lbl_numfp->setWordWrap(true);
        spinBox_numGPUs = new QSpinBox(centralWidget);
        spinBox_numGPUs->setObjectName(QStringLiteral("spinBox_numGPUs"));
        spinBox_numGPUs->setGeometry(QRect(343, 208, 48, 24));
        spinBox_numGPUs->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lbl_numgpus = new QLabel(centralWidget);
        lbl_numgpus->setObjectName(QStringLiteral("lbl_numgpus"));
        lbl_numgpus->setGeometry(QRect(269, 208, 66, 16));
        lineEdit_potbound = new QLineEdit(centralWidget);
        lineEdit_potbound->setObjectName(QStringLiteral("lineEdit_potbound"));
        lineEdit_potbound->setGeometry(QRect(121, 175, 125, 21));
        sizePolicy1.setHeightForWidth(lineEdit_potbound->sizePolicy().hasHeightForWidth());
        lineEdit_potbound->setSizePolicy(sizePolicy1);
        lineEdit_potbound->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lbl_potBound = new QLabel(centralWidget);
        lbl_potBound->setObjectName(QStringLiteral("lbl_potBound"));
        lbl_potBound->setGeometry(QRect(17, 175, 96, 16));
        label_sliceThickness = new QLabel(centralWidget);
        label_sliceThickness->setObjectName(QStringLiteral("label_sliceThickness"));
        label_sliceThickness->setGeometry(QRect(17, 244, 94, 16));
        lineEdit_sliceThickness = new QLineEdit(centralWidget);
        lineEdit_sliceThickness->setObjectName(QStringLiteral("lineEdit_sliceThickness"));
        lineEdit_sliceThickness->setGeometry(QRect(119, 244, 125, 21));
        sizePolicy1.setHeightForWidth(lineEdit_sliceThickness->sizePolicy().hasHeightForWidth());
        lineEdit_sliceThickness->setSizePolicy(sizePolicy1);
        lineEdit_sliceThickness->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        spinBox_numThreads = new QSpinBox(centralWidget);
        spinBox_numThreads->setObjectName(QStringLiteral("spinBox_numThreads"));
        spinBox_numThreads->setGeometry(QRect(359, 244, 48, 24));
        spinBox_numThreads->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lbl_numthreads = new QLabel(centralWidget);
        lbl_numthreads->setObjectName(QStringLiteral("lbl_numthreads"));
        lbl_numthreads->setGeometry(QRect(269, 244, 82, 16));
        lbl_E0 = new QLabel(centralWidget);
        lbl_E0->setObjectName(QStringLiteral("lbl_E0"));
        lbl_E0->setGeometry(QRect(269, 175, 77, 16));
        lineEdit_E0 = new QLineEdit(centralWidget);
        lineEdit_E0->setObjectName(QStringLiteral("lineEdit_E0"));
        lineEdit_E0->setGeometry(QRect(354, 175, 125, 21));
        sizePolicy1.setHeightForWidth(lineEdit_E0->sizePolicy().hasHeightForWidth());
        lineEdit_E0->setSizePolicy(sizePolicy1);
        lineEdit_E0->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        PRISMMainWindow->setCentralWidget(centralWidget);
        box_simulationSettings->raise();
        box_calculationSettings->raise();
        box_sampleSettings->raise();
        box_sampleSettings->raise();
        menuBar = new QMenuBar(PRISMMainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 901, 22));
        PRISMMainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(PRISMMainWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        PRISMMainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(PRISMMainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        PRISMMainWindow->setStatusBar(statusBar);

        retranslateUi(PRISMMainWindow);

        QMetaObject::connectSlotsByName(PRISMMainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *PRISMMainWindow)
    {
        PRISMMainWindow->setWindowTitle(QApplication::translate("PRISMMainWindow", "PRISMMainWindow", Q_NULLPTR));
        box_simulationSettings->setTitle(QApplication::translate("PRISMMainWindow", "Simulation Settings", Q_NULLPTR));
        box_calculationSettings->setTitle(QApplication::translate("PRISMMainWindow", "Calculation Settings", Q_NULLPTR));
        box_sampleSettings->setTitle(QApplication::translate("PRISMMainWindow", "Sample Settings", Q_NULLPTR));
        btn_loadCoordinates->setText(QApplication::translate("PRISMMainWindow", "Load Coords", Q_NULLPTR));
        btn_saveCoordinates->setText(QApplication::translate("PRISMMainWindow", "Save Coords", Q_NULLPTR));
        label_8->setText(QString());
        label_cellDim->setText(QApplication::translate("PRISMMainWindow", "Cell Dim.", Q_NULLPTR));
        label_6->setText(QApplication::translate("PRISMMainWindow", "Tile Cells", Q_NULLPTR));
        label_2->setText(QApplication::translate("PRISMMainWindow", "X", Q_NULLPTR));
        lineEdit_cellDimX->setText(QApplication::translate("PRISMMainWindow", "100", Q_NULLPTR));
        label_3->setText(QApplication::translate("PRISMMainWindow", "Y", Q_NULLPTR));
        lineEdit_cellDimY->setText(QApplication::translate("PRISMMainWindow", "100", Q_NULLPTR));
        label_4->setText(QApplication::translate("PRISMMainWindow", "Z", Q_NULLPTR));
        lineEdit_cellDimZ->setText(QApplication::translate("PRISMMainWindow", "80", Q_NULLPTR));
        label_9->setText(QString());
        label_5->setText(QApplication::translate("PRISMMainWindow", "A", Q_NULLPTR));
        label_7->setText(QApplication::translate("PRISMMainWindow", "UCs", Q_NULLPTR));
        lbl_algo->setText(QApplication::translate("PRISMMainWindow", "Algorithm", Q_NULLPTR));
        radBtn_PRISM->setText(QApplication::translate("PRISMMainWindow", "PRISM", Q_NULLPTR));
        radBtn_Multislice->setText(QApplication::translate("PRISMMainWindow", "Multislice", Q_NULLPTR));
        btn_go->setText(QApplication::translate("PRISMMainWindow", "Calculate!", Q_NULLPTR));
        btn_outputfile_browse->setText(QApplication::translate("PRISMMainWindow", "Browse", Q_NULLPTR));
        lbl_outputfile->setText(QApplication::translate("PRISMMainWindow", "Output File", Q_NULLPTR));
        btn_atomsfile_browse->setText(QApplication::translate("PRISMMainWindow", "Browse", Q_NULLPTR));
        lbl_atomsfile->setText(QApplication::translate("PRISMMainWindow", "Atoms File ", Q_NULLPTR));
        label->setText(QApplication::translate("PRISMMainWindow", "Interpolation factor (f)", Q_NULLPTR));
        lbl_pixelSize->setText(QApplication::translate("PRISMMainWindow", "Pixel Size", Q_NULLPTR));
        lineedit_pixelSize->setText(QApplication::translate("PRISMMainWindow", "0.1", Q_NULLPTR));
        lineEdit_alphaBeamMax->setText(QApplication::translate("PRISMMainWindow", "24", Q_NULLPTR));
        label_alphaBeamMax->setText(QApplication::translate("PRISMMainWindow", "Max Angle (mrad)", Q_NULLPTR));
        lbl_numfp->setText(QApplication::translate("PRISMMainWindow", "Number of Frozen Phonons", Q_NULLPTR));
        lbl_numgpus->setText(QApplication::translate("PRISMMainWindow", "Num GPUs", Q_NULLPTR));
        lineEdit_potbound->setText(QApplication::translate("PRISMMainWindow", "1", Q_NULLPTR));
        lbl_potBound->setText(QApplication::translate("PRISMMainWindow", "Potential Bound", Q_NULLPTR));
        label_sliceThickness->setText(QApplication::translate("PRISMMainWindow", "Slice Thickness", Q_NULLPTR));
        lineEdit_sliceThickness->setText(QApplication::translate("PRISMMainWindow", "2", Q_NULLPTR));
        lbl_numthreads->setText(QApplication::translate("PRISMMainWindow", "Num Threads", Q_NULLPTR));
        lbl_E0->setText(QApplication::translate("PRISMMainWindow", "Energy (keV)", Q_NULLPTR));
        lineEdit_E0->setText(QApplication::translate("PRISMMainWindow", "80", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class PRISMMainWindow: public Ui_PRISMMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PRISMMAINWINDOW_H
