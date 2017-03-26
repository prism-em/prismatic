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
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
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
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout_16;
    QGroupBox *box_sampleSettings;
    QHBoxLayout *horizontalLayout_13;
    QVBoxLayout *verticalLayout_13;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *btn_atomsfile_browse;
    QPushButton *btn_saveCoordinates;
    QHBoxLayout *horizontalLayout_12;
    QVBoxLayout *verticalLayout_12;
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
    QGroupBox *box_simulationSettings;
    QHBoxLayout *horizontalLayout_11;
    QVBoxLayout *verticalLayout_10;
    QHBoxLayout *horizontalLayout_4;
    QVBoxLayout *verticalLayout_5;
    QHBoxLayout *horizontalLayout_5;
    QLabel *lbl_pixelSize;
    QLineEdit *lineedit_pixelSize;
    QHBoxLayout *horizontalLayout_3;
    QLabel *lbl_E0;
    QLineEdit *lineEdit_E0;
    QHBoxLayout *horizontalLayout_7;
    QLabel *lbl_potBound;
    QLineEdit *lineEdit_potbound;
    QHBoxLayout *horizontalLayout;
    QLabel *label_alphaBeamMax;
    QLineEdit *lineEdit_alphaBeamMax;
    QVBoxLayout *verticalLayout_6;
    QLabel *label_10;
    QLabel *label_11;
    QHBoxLayout *horizontalLayout_6;
    QLabel *lbl_numfp;
    QSpinBox *spinBox_numFP;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_sliceThickness;
    QLineEdit *lineEdit_sliceThickness;
    QHBoxLayout *horizontalLayout_10;
    QVBoxLayout *verticalLayout_9;
    QLabel *label;
    QHBoxLayout *horizontalLayout_8;
    QVBoxLayout *verticalLayout_7;
    QLabel *label_13;
    QLineEdit *lineedit_interpFactor_x;
    QVBoxLayout *verticalLayout_8;
    QLabel *label_14;
    QLineEdit *lineedit_interpFactor_y;
    QGroupBox *box_calculationSettings;
    QHBoxLayout *horizontalLayout_20;
    QVBoxLayout *verticalLayout_14;
    QHBoxLayout *horizontalLayout_14;
    QLabel *lbl_outputfile;
    QLineEdit *lineedit_outputfile;
    QHBoxLayout *horizontalLayout_19;
    QCheckBox *checkBox;
    QLabel *label_12;
    QHBoxLayout *horizontalLayout_18;
    QRadioButton *radioButton;
    QLineEdit *lineEdit;
    QSpacerItem *horizontalSpacer;
    QHBoxLayout *horizontalLayout_17;
    QRadioButton *radioButton_2;
    QLineEdit *lineEdit_2;
    QSpacerItem *horizontalSpacer_2;
    QHBoxLayout *horizontalLayout_16;
    QLabel *lbl_numthreads;
    QSpinBox *spinBox_numThreads;
    QHBoxLayout *horizontalLayout_15;
    QLabel *lbl_numgpus;
    QSpinBox *spinBox_numGPUs;
    QSpacerItem *horizontalSpacer_3;
    QHBoxLayout *horizontalLayout_22;
    QHBoxLayout *horizontalLayout_21;
    QLabel *lbl_algo;
    QRadioButton *radBtn_PRISM;
    QRadioButton *radBtn_Multislice;
    QSpacerItem *horizontalSpacer_4;
    QVBoxLayout *verticalLayout_15;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QPushButton *btn_go;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PRISMMainWindow)
    {
        if (PRISMMainWindow->objectName().isEmpty())
            PRISMMainWindow->setObjectName(QStringLiteral("PRISMMainWindow"));
        PRISMMainWindow->resize(333, 849);
        centralWidget = new QWidget(PRISMMainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        verticalLayout_16 = new QVBoxLayout();
        verticalLayout_16->setSpacing(6);
        verticalLayout_16->setObjectName(QStringLiteral("verticalLayout_16"));
        box_sampleSettings = new QGroupBox(centralWidget);
        box_sampleSettings->setObjectName(QStringLiteral("box_sampleSettings"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(box_sampleSettings->sizePolicy().hasHeightForWidth());
        box_sampleSettings->setSizePolicy(sizePolicy);
        horizontalLayout_13 = new QHBoxLayout(box_sampleSettings);
        horizontalLayout_13->setSpacing(0);
        horizontalLayout_13->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_13->setObjectName(QStringLiteral("horizontalLayout_13"));
        horizontalLayout_13->setContentsMargins(2, 10, 2, 2);
        verticalLayout_13 = new QVBoxLayout();
        verticalLayout_13->setSpacing(0);
        verticalLayout_13->setObjectName(QStringLiteral("verticalLayout_13"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(0);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        btn_atomsfile_browse = new QPushButton(box_sampleSettings);
        btn_atomsfile_browse->setObjectName(QStringLiteral("btn_atomsfile_browse"));

        horizontalLayout_2->addWidget(btn_atomsfile_browse);

        btn_saveCoordinates = new QPushButton(box_sampleSettings);
        btn_saveCoordinates->setObjectName(QStringLiteral("btn_saveCoordinates"));
        QSizePolicy sizePolicy1(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(btn_saveCoordinates->sizePolicy().hasHeightForWidth());
        btn_saveCoordinates->setSizePolicy(sizePolicy1);

        horizontalLayout_2->addWidget(btn_saveCoordinates);


        verticalLayout_13->addLayout(horizontalLayout_2);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setSpacing(6);
        horizontalLayout_12->setObjectName(QStringLiteral("horizontalLayout_12"));
        verticalLayout_12 = new QVBoxLayout();
        verticalLayout_12->setSpacing(6);
        verticalLayout_12->setObjectName(QStringLiteral("verticalLayout_12"));
        label_8 = new QLabel(box_sampleSettings);
        label_8->setObjectName(QStringLiteral("label_8"));

        verticalLayout_12->addWidget(label_8);

        label_cellDim = new QLabel(box_sampleSettings);
        label_cellDim->setObjectName(QStringLiteral("label_cellDim"));
        label_cellDim->setAlignment(Qt::AlignCenter);

        verticalLayout_12->addWidget(label_cellDim);

        label_6 = new QLabel(box_sampleSettings);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setAlignment(Qt::AlignCenter);

        verticalLayout_12->addWidget(label_6);


        horizontalLayout_12->addLayout(verticalLayout_12);

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


        horizontalLayout_12->addLayout(verticalLayout);

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


        horizontalLayout_12->addLayout(verticalLayout_2);

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


        horizontalLayout_12->addLayout(verticalLayout_3);

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


        horizontalLayout_12->addLayout(verticalLayout_4);


        verticalLayout_13->addLayout(horizontalLayout_12);


        horizontalLayout_13->addLayout(verticalLayout_13);


        verticalLayout_16->addWidget(box_sampleSettings);

        box_simulationSettings = new QGroupBox(centralWidget);
        box_simulationSettings->setObjectName(QStringLiteral("box_simulationSettings"));
        horizontalLayout_11 = new QHBoxLayout(box_simulationSettings);
        horizontalLayout_11->setSpacing(0);
        horizontalLayout_11->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_11->setObjectName(QStringLiteral("horizontalLayout_11"));
        horizontalLayout_11->setContentsMargins(2, 10, 2, 2);
        verticalLayout_10 = new QVBoxLayout();
        verticalLayout_10->setSpacing(0);
        verticalLayout_10->setObjectName(QStringLiteral("verticalLayout_10"));
        verticalLayout_10->setContentsMargins(-1, 10, -1, -1);
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(30);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, -1, -1, -1);
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(0);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        lbl_pixelSize = new QLabel(box_simulationSettings);
        lbl_pixelSize->setObjectName(QStringLiteral("lbl_pixelSize"));

        horizontalLayout_5->addWidget(lbl_pixelSize);

        lineedit_pixelSize = new QLineEdit(box_simulationSettings);
        lineedit_pixelSize->setObjectName(QStringLiteral("lineedit_pixelSize"));
        sizePolicy1.setHeightForWidth(lineedit_pixelSize->sizePolicy().hasHeightForWidth());
        lineedit_pixelSize->setSizePolicy(sizePolicy1);
        lineedit_pixelSize->setMaximumSize(QSize(30, 16777215));
        lineedit_pixelSize->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_5->addWidget(lineedit_pixelSize);


        verticalLayout_5->addLayout(horizontalLayout_5);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        lbl_E0 = new QLabel(box_simulationSettings);
        lbl_E0->setObjectName(QStringLiteral("lbl_E0"));

        horizontalLayout_3->addWidget(lbl_E0);

        lineEdit_E0 = new QLineEdit(box_simulationSettings);
        lineEdit_E0->setObjectName(QStringLiteral("lineEdit_E0"));
        sizePolicy1.setHeightForWidth(lineEdit_E0->sizePolicy().hasHeightForWidth());
        lineEdit_E0->setSizePolicy(sizePolicy1);
        lineEdit_E0->setMaximumSize(QSize(30, 16777215));
        lineEdit_E0->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_3->addWidget(lineEdit_E0);


        verticalLayout_5->addLayout(horizontalLayout_3);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        lbl_potBound = new QLabel(box_simulationSettings);
        lbl_potBound->setObjectName(QStringLiteral("lbl_potBound"));

        horizontalLayout_7->addWidget(lbl_potBound);

        lineEdit_potbound = new QLineEdit(box_simulationSettings);
        lineEdit_potbound->setObjectName(QStringLiteral("lineEdit_potbound"));
        sizePolicy1.setHeightForWidth(lineEdit_potbound->sizePolicy().hasHeightForWidth());
        lineEdit_potbound->setSizePolicy(sizePolicy1);
        lineEdit_potbound->setMaximumSize(QSize(30, 16777215));
        lineEdit_potbound->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_7->addWidget(lineEdit_potbound);


        verticalLayout_5->addLayout(horizontalLayout_7);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_alphaBeamMax = new QLabel(box_simulationSettings);
        label_alphaBeamMax->setObjectName(QStringLiteral("label_alphaBeamMax"));
        label_alphaBeamMax->setMinimumSize(QSize(0, 40));
        label_alphaBeamMax->setWordWrap(true);

        horizontalLayout->addWidget(label_alphaBeamMax);

        lineEdit_alphaBeamMax = new QLineEdit(box_simulationSettings);
        lineEdit_alphaBeamMax->setObjectName(QStringLiteral("lineEdit_alphaBeamMax"));
        sizePolicy1.setHeightForWidth(lineEdit_alphaBeamMax->sizePolicy().hasHeightForWidth());
        lineEdit_alphaBeamMax->setSizePolicy(sizePolicy1);
        lineEdit_alphaBeamMax->setMaximumSize(QSize(30, 16777215));
        lineEdit_alphaBeamMax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout->addWidget(lineEdit_alphaBeamMax);


        verticalLayout_5->addLayout(horizontalLayout);


        horizontalLayout_4->addLayout(verticalLayout_5);

        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        label_10 = new QLabel(box_simulationSettings);
        label_10->setObjectName(QStringLiteral("label_10"));
        label_10->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_6->addWidget(label_10);

        label_11 = new QLabel(box_simulationSettings);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_6->addWidget(label_11);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        lbl_numfp = new QLabel(box_simulationSettings);
        lbl_numfp->setObjectName(QStringLiteral("lbl_numfp"));
        lbl_numfp->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lbl_numfp->setWordWrap(true);

        horizontalLayout_6->addWidget(lbl_numfp);

        spinBox_numFP = new QSpinBox(box_simulationSettings);
        spinBox_numFP->setObjectName(QStringLiteral("spinBox_numFP"));
        spinBox_numFP->setMaximumSize(QSize(40, 16777215));
        spinBox_numFP->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_6->addWidget(spinBox_numFP);


        verticalLayout_6->addLayout(horizontalLayout_6);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        label_sliceThickness = new QLabel(box_simulationSettings);
        label_sliceThickness->setObjectName(QStringLiteral("label_sliceThickness"));
        label_sliceThickness->setMinimumSize(QSize(0, 40));
        label_sliceThickness->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_sliceThickness->setWordWrap(true);

        horizontalLayout_9->addWidget(label_sliceThickness);

        lineEdit_sliceThickness = new QLineEdit(box_simulationSettings);
        lineEdit_sliceThickness->setObjectName(QStringLiteral("lineEdit_sliceThickness"));
        sizePolicy1.setHeightForWidth(lineEdit_sliceThickness->sizePolicy().hasHeightForWidth());
        lineEdit_sliceThickness->setSizePolicy(sizePolicy1);
        lineEdit_sliceThickness->setMaximumSize(QSize(30, 16777215));
        lineEdit_sliceThickness->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_9->addWidget(lineEdit_sliceThickness);


        verticalLayout_6->addLayout(horizontalLayout_9);


        horizontalLayout_4->addLayout(verticalLayout_6);


        verticalLayout_10->addLayout(horizontalLayout_4);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(0);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        verticalLayout_9 = new QVBoxLayout();
        verticalLayout_9->setSpacing(0);
        verticalLayout_9->setObjectName(QStringLiteral("verticalLayout_9"));
        verticalLayout_9->setContentsMargins(-1, 20, -1, -1);
        label = new QLabel(box_simulationSettings);
        label->setObjectName(QStringLiteral("label"));
        label->setMaximumSize(QSize(16777215, 30));
        label->setWordWrap(true);

        verticalLayout_9->addWidget(label);


        horizontalLayout_10->addLayout(verticalLayout_9);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        verticalLayout_7 = new QVBoxLayout();
        verticalLayout_7->setSpacing(0);
        verticalLayout_7->setObjectName(QStringLiteral("verticalLayout_7"));
        label_13 = new QLabel(box_simulationSettings);
        label_13->setObjectName(QStringLiteral("label_13"));
        label_13->setAlignment(Qt::AlignCenter);

        verticalLayout_7->addWidget(label_13);

        lineedit_interpFactor_x = new QLineEdit(box_simulationSettings);
        lineedit_interpFactor_x->setObjectName(QStringLiteral("lineedit_interpFactor_x"));
        lineedit_interpFactor_x->setMaximumSize(QSize(30, 16777215));

        verticalLayout_7->addWidget(lineedit_interpFactor_x);


        horizontalLayout_8->addLayout(verticalLayout_7);

        verticalLayout_8 = new QVBoxLayout();
        verticalLayout_8->setSpacing(0);
        verticalLayout_8->setObjectName(QStringLiteral("verticalLayout_8"));
        label_14 = new QLabel(box_simulationSettings);
        label_14->setObjectName(QStringLiteral("label_14"));
        label_14->setAlignment(Qt::AlignCenter);

        verticalLayout_8->addWidget(label_14);

        lineedit_interpFactor_y = new QLineEdit(box_simulationSettings);
        lineedit_interpFactor_y->setObjectName(QStringLiteral("lineedit_interpFactor_y"));
        lineedit_interpFactor_y->setMaximumSize(QSize(30, 16777215));

        verticalLayout_8->addWidget(lineedit_interpFactor_y);


        horizontalLayout_8->addLayout(verticalLayout_8);


        horizontalLayout_10->addLayout(horizontalLayout_8);


        verticalLayout_10->addLayout(horizontalLayout_10);


        horizontalLayout_11->addLayout(verticalLayout_10);


        verticalLayout_16->addWidget(box_simulationSettings);

        box_calculationSettings = new QGroupBox(centralWidget);
        box_calculationSettings->setObjectName(QStringLiteral("box_calculationSettings"));
        horizontalLayout_20 = new QHBoxLayout(box_calculationSettings);
        horizontalLayout_20->setSpacing(0);
        horizontalLayout_20->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_20->setObjectName(QStringLiteral("horizontalLayout_20"));
        horizontalLayout_20->setContentsMargins(2, 10, 2, 2);
        verticalLayout_14 = new QVBoxLayout();
        verticalLayout_14->setSpacing(6);
        verticalLayout_14->setObjectName(QStringLiteral("verticalLayout_14"));
        verticalLayout_14->setContentsMargins(-1, 10, -1, -1);
        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setSpacing(6);
        horizontalLayout_14->setObjectName(QStringLiteral("horizontalLayout_14"));
        lbl_outputfile = new QLabel(box_calculationSettings);
        lbl_outputfile->setObjectName(QStringLiteral("lbl_outputfile"));

        horizontalLayout_14->addWidget(lbl_outputfile);

        lineedit_outputfile = new QLineEdit(box_calculationSettings);
        lineedit_outputfile->setObjectName(QStringLiteral("lineedit_outputfile"));

        horizontalLayout_14->addWidget(lineedit_outputfile);


        verticalLayout_14->addLayout(horizontalLayout_14);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setSpacing(6);
        horizontalLayout_19->setObjectName(QStringLiteral("horizontalLayout_19"));
        checkBox = new QCheckBox(box_calculationSettings);
        checkBox->setObjectName(QStringLiteral("checkBox"));

        horizontalLayout_19->addWidget(checkBox);

        label_12 = new QLabel(box_calculationSettings);
        label_12->setObjectName(QStringLiteral("label_12"));

        horizontalLayout_19->addWidget(label_12);


        verticalLayout_14->addLayout(horizontalLayout_19);

        horizontalLayout_18 = new QHBoxLayout();
        horizontalLayout_18->setSpacing(6);
        horizontalLayout_18->setObjectName(QStringLiteral("horizontalLayout_18"));
        radioButton = new QRadioButton(box_calculationSettings);
        radioButton->setObjectName(QStringLiteral("radioButton"));

        horizontalLayout_18->addWidget(radioButton);

        lineEdit = new QLineEdit(box_calculationSettings);
        lineEdit->setObjectName(QStringLiteral("lineEdit"));
        lineEdit->setMaximumSize(QSize(30, 16777215));

        horizontalLayout_18->addWidget(lineEdit);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_18->addItem(horizontalSpacer);


        verticalLayout_14->addLayout(horizontalLayout_18);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setSpacing(6);
        horizontalLayout_17->setObjectName(QStringLiteral("horizontalLayout_17"));
        radioButton_2 = new QRadioButton(box_calculationSettings);
        radioButton_2->setObjectName(QStringLiteral("radioButton_2"));

        horizontalLayout_17->addWidget(radioButton_2);

        lineEdit_2 = new QLineEdit(box_calculationSettings);
        lineEdit_2->setObjectName(QStringLiteral("lineEdit_2"));
        lineEdit_2->setMaximumSize(QSize(30, 16777215));

        horizontalLayout_17->addWidget(lineEdit_2);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_17->addItem(horizontalSpacer_2);


        verticalLayout_14->addLayout(horizontalLayout_17);

        horizontalLayout_16 = new QHBoxLayout();
        horizontalLayout_16->setSpacing(6);
        horizontalLayout_16->setObjectName(QStringLiteral("horizontalLayout_16"));
        lbl_numthreads = new QLabel(box_calculationSettings);
        lbl_numthreads->setObjectName(QStringLiteral("lbl_numthreads"));

        horizontalLayout_16->addWidget(lbl_numthreads);

        spinBox_numThreads = new QSpinBox(box_calculationSettings);
        spinBox_numThreads->setObjectName(QStringLiteral("spinBox_numThreads"));
        spinBox_numThreads->setMaximumSize(QSize(40, 16777215));
        spinBox_numThreads->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_16->addWidget(spinBox_numThreads);

        horizontalLayout_15 = new QHBoxLayout();
        horizontalLayout_15->setSpacing(6);
        horizontalLayout_15->setObjectName(QStringLiteral("horizontalLayout_15"));
        lbl_numgpus = new QLabel(box_calculationSettings);
        lbl_numgpus->setObjectName(QStringLiteral("lbl_numgpus"));

        horizontalLayout_15->addWidget(lbl_numgpus);

        spinBox_numGPUs = new QSpinBox(box_calculationSettings);
        spinBox_numGPUs->setObjectName(QStringLiteral("spinBox_numGPUs"));
        spinBox_numGPUs->setMaximumSize(QSize(40, 16777215));
        spinBox_numGPUs->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_15->addWidget(spinBox_numGPUs);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_15->addItem(horizontalSpacer_3);


        horizontalLayout_16->addLayout(horizontalLayout_15);


        verticalLayout_14->addLayout(horizontalLayout_16);

        horizontalLayout_22 = new QHBoxLayout();
        horizontalLayout_22->setSpacing(6);
        horizontalLayout_22->setObjectName(QStringLiteral("horizontalLayout_22"));
        horizontalLayout_21 = new QHBoxLayout();
        horizontalLayout_21->setSpacing(6);
        horizontalLayout_21->setObjectName(QStringLiteral("horizontalLayout_21"));
        lbl_algo = new QLabel(box_calculationSettings);
        lbl_algo->setObjectName(QStringLiteral("lbl_algo"));
        lbl_algo->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_21->addWidget(lbl_algo);

        radBtn_PRISM = new QRadioButton(box_calculationSettings);
        radBtn_PRISM->setObjectName(QStringLiteral("radBtn_PRISM"));

        horizontalLayout_21->addWidget(radBtn_PRISM);

        radBtn_Multislice = new QRadioButton(box_calculationSettings);
        radBtn_Multislice->setObjectName(QStringLiteral("radBtn_Multislice"));

        horizontalLayout_21->addWidget(radBtn_Multislice);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_21->addItem(horizontalSpacer_4);


        horizontalLayout_22->addLayout(horizontalLayout_21);


        verticalLayout_14->addLayout(horizontalLayout_22);

        verticalLayout_15 = new QVBoxLayout();
        verticalLayout_15->setSpacing(6);
        verticalLayout_15->setObjectName(QStringLiteral("verticalLayout_15"));
        pushButton = new QPushButton(box_calculationSettings);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setMinimumSize(QSize(0, 50));
        QFont font;
        font.setPointSize(14);
        pushButton->setFont(font);

        verticalLayout_15->addWidget(pushButton);

        pushButton_2 = new QPushButton(box_calculationSettings);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));
        pushButton_2->setMinimumSize(QSize(0, 50));
        pushButton_2->setFont(font);

        verticalLayout_15->addWidget(pushButton_2);

        btn_go = new QPushButton(box_calculationSettings);
        btn_go->setObjectName(QStringLiteral("btn_go"));
        btn_go->setMinimumSize(QSize(0, 50));
        btn_go->setFont(font);

        verticalLayout_15->addWidget(btn_go);


        verticalLayout_14->addLayout(verticalLayout_15);


        horizontalLayout_20->addLayout(verticalLayout_14);


        verticalLayout_16->addWidget(box_calculationSettings);


        gridLayout->addLayout(verticalLayout_16, 0, 0, 1, 1);

        PRISMMainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(PRISMMainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 333, 22));
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
        box_sampleSettings->setTitle(QApplication::translate("PRISMMainWindow", "Sample Settings", Q_NULLPTR));
        btn_atomsfile_browse->setText(QApplication::translate("PRISMMainWindow", "Load Coords", Q_NULLPTR));
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
        box_simulationSettings->setTitle(QApplication::translate("PRISMMainWindow", "Simulation Settings", Q_NULLPTR));
        lbl_pixelSize->setText(QApplication::translate("PRISMMainWindow", "Pixel Size (A)", Q_NULLPTR));
        lineedit_pixelSize->setText(QApplication::translate("PRISMMainWindow", "0.1", Q_NULLPTR));
        lbl_E0->setText(QApplication::translate("PRISMMainWindow", "Energy (keV)", Q_NULLPTR));
        lineEdit_E0->setText(QApplication::translate("PRISMMainWindow", "80", Q_NULLPTR));
        lbl_potBound->setText(QApplication::translate("PRISMMainWindow", "Potential Bound (A)", Q_NULLPTR));
        lineEdit_potbound->setText(QApplication::translate("PRISMMainWindow", "1", Q_NULLPTR));
        label_alphaBeamMax->setText(QApplication::translate("PRISMMainWindow", "Probe Semiangle (mrads)", Q_NULLPTR));
        lineEdit_alphaBeamMax->setText(QApplication::translate("PRISMMainWindow", "20", Q_NULLPTR));
        label_10->setText(QApplication::translate("PRISMMainWindow", "l = 0.1 A", Q_NULLPTR));
        label_11->setText(QApplication::translate("PRISMMainWindow", "a_max = 24", Q_NULLPTR));
        lbl_numfp->setText(QApplication::translate("PRISMMainWindow", "# of FP", Q_NULLPTR));
        label_sliceThickness->setText(QApplication::translate("PRISMMainWindow", "Slice Thickness (A)", Q_NULLPTR));
        lineEdit_sliceThickness->setText(QApplication::translate("PRISMMainWindow", "2", Q_NULLPTR));
        label->setText(QApplication::translate("PRISMMainWindow", "PRISM Interpolation factors", Q_NULLPTR));
        label_13->setText(QApplication::translate("PRISMMainWindow", "X", Q_NULLPTR));
        label_14->setText(QApplication::translate("PRISMMainWindow", "Y", Q_NULLPTR));
        box_calculationSettings->setTitle(QApplication::translate("PRISMMainWindow", "Calculation Settings", Q_NULLPTR));
        lbl_outputfile->setText(QApplication::translate("PRISMMainWindow", "Output File", Q_NULLPTR));
        checkBox->setText(QApplication::translate("PRISMMainWindow", "Save projected potentials", Q_NULLPTR));
        label_12->setText(QString());
        radioButton->setText(QApplication::translate("PRISMMainWindow", "3D output, step size (mrad)", Q_NULLPTR));
        radioButton_2->setText(QApplication::translate("PRISMMainWindow", "4D output, step size (mrad)", Q_NULLPTR));
        lbl_numthreads->setText(QApplication::translate("PRISMMainWindow", "Num Threads", Q_NULLPTR));
        lbl_numgpus->setText(QApplication::translate("PRISMMainWindow", "Num GPUs", Q_NULLPTR));
        lbl_algo->setText(QApplication::translate("PRISMMainWindow", "Algorithm", Q_NULLPTR));
        radBtn_PRISM->setText(QApplication::translate("PRISMMainWindow", "PRISM", Q_NULLPTR));
        radBtn_Multislice->setText(QApplication::translate("PRISMMainWindow", "Multislice", Q_NULLPTR));
        pushButton->setText(QApplication::translate("PRISMMainWindow", "Calculate Potentials", Q_NULLPTR));
        pushButton_2->setText(QApplication::translate("PRISMMainWindow", "Calculate S-matrix", Q_NULLPTR));
        btn_go->setText(QApplication::translate("PRISMMainWindow", "Full Calculation", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class PRISMMainWindow: public Ui_PRISMMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PRISMMAINWINDOW_H
