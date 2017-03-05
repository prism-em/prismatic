/********************************************************************************
** Form generated from reading UI file 'prismmainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PRISMMAINWINDOW_H
#define UI_PRISMMAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
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
    QVBoxLayout *verticalLayout_3;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *lbl_atomsfile;
    QLineEdit *lineedit_atomsfile;
    QPushButton *btn_atomsfile_browse;
    QHBoxLayout *horizontalLayout_3;
    QLabel *lbl_outputfile;
    QLineEdit *lineedit_outputfile;
    QPushButton *btn_outputfile_browse;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QLineEdit *lineedit_f;
    QSpacerItem *horizontalSpacer;
    QHBoxLayout *horizontalLayout_4;
    QVBoxLayout *verticalLayout_2;
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout_5;
    QLabel *lbl_pixelSize;
    QLineEdit *lineedit_pixelSize;
    QHBoxLayout *horizontalLayout_7;
    QLabel *lbl_numfp;
    QSpinBox *spinBox_numFP;
    QHBoxLayout *horizontalLayout_6;
    QLabel *lbl_potBound;
    QLineEdit *lineEdit_potbound;
    QHBoxLayout *horizontalLayout_10;
    QLabel *lbl_E0;
    QLineEdit *lineEdit_E0;
    QHBoxLayout *horizontalLayout_12;
    QLabel *label_alphaBeamMax;
    QLineEdit *lineEdit_alphaBeamMax;
    QHBoxLayout *horizontalLayout_9;
    QLabel *lbl_numgpus;
    QSpinBox *spinBox_numGPUs;
    QHBoxLayout *horizontalLayout_13;
    QLabel *label_sliceThickness;
    QLineEdit *lineEdit_sliceThickness;
    QHBoxLayout *horizontalLayout_8;
    QLabel *lbl_numthreads;
    QSpinBox *spinBox_numThreads;
    QHBoxLayout *horizontalLayout_11;
    QLabel *label_cellDim;
    QLineEdit *lineEdit_cellDimX;
    QLineEdit *lineEdit_cellDimY;
    QLineEdit *lineEdit_cellDimZ;
    QPushButton *btn_go;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PRISMMainWindow)
    {
        if (PRISMMainWindow->objectName().isEmpty())
            PRISMMainWindow->setObjectName(QStringLiteral("PRISMMainWindow"));
        PRISMMainWindow->resize(778, 347);
        centralWidget = new QWidget(PRISMMainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        verticalLayout_3 = new QVBoxLayout(centralWidget);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        lbl_atomsfile = new QLabel(centralWidget);
        lbl_atomsfile->setObjectName(QStringLiteral("lbl_atomsfile"));

        horizontalLayout->addWidget(lbl_atomsfile);

        lineedit_atomsfile = new QLineEdit(centralWidget);
        lineedit_atomsfile->setObjectName(QStringLiteral("lineedit_atomsfile"));

        horizontalLayout->addWidget(lineedit_atomsfile);

        btn_atomsfile_browse = new QPushButton(centralWidget);
        btn_atomsfile_browse->setObjectName(QStringLiteral("btn_atomsfile_browse"));

        horizontalLayout->addWidget(btn_atomsfile_browse);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        lbl_outputfile = new QLabel(centralWidget);
        lbl_outputfile->setObjectName(QStringLiteral("lbl_outputfile"));

        horizontalLayout_3->addWidget(lbl_outputfile);

        lineedit_outputfile = new QLineEdit(centralWidget);
        lineedit_outputfile->setObjectName(QStringLiteral("lineedit_outputfile"));

        horizontalLayout_3->addWidget(lineedit_outputfile);

        btn_outputfile_browse = new QPushButton(centralWidget);
        btn_outputfile_browse->setObjectName(QStringLiteral("btn_outputfile_browse"));

        horizontalLayout_3->addWidget(btn_outputfile_browse);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout_2->addWidget(label);

        lineedit_f = new QLineEdit(centralWidget);
        lineedit_f->setObjectName(QStringLiteral("lineedit_f"));
        lineedit_f->setMaximumSize(QSize(50, 16777215));

        horizontalLayout_2->addWidget(lineedit_f);

        horizontalSpacer = new QSpacerItem(400, 20, QSizePolicy::Maximum, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        lbl_pixelSize = new QLabel(centralWidget);
        lbl_pixelSize->setObjectName(QStringLiteral("lbl_pixelSize"));

        horizontalLayout_5->addWidget(lbl_pixelSize);

        lineedit_pixelSize = new QLineEdit(centralWidget);
        lineedit_pixelSize->setObjectName(QStringLiteral("lineedit_pixelSize"));
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(lineedit_pixelSize->sizePolicy().hasHeightForWidth());
        lineedit_pixelSize->setSizePolicy(sizePolicy);
        lineedit_pixelSize->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_5->addWidget(lineedit_pixelSize);


        gridLayout->addLayout(horizontalLayout_5, 0, 0, 1, 1);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        lbl_numfp = new QLabel(centralWidget);
        lbl_numfp->setObjectName(QStringLiteral("lbl_numfp"));
        lbl_numfp->setWordWrap(true);

        horizontalLayout_7->addWidget(lbl_numfp);

        spinBox_numFP = new QSpinBox(centralWidget);
        spinBox_numFP->setObjectName(QStringLiteral("spinBox_numFP"));
        spinBox_numFP->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_7->addWidget(spinBox_numFP);


        gridLayout->addLayout(horizontalLayout_7, 0, 1, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        lbl_potBound = new QLabel(centralWidget);
        lbl_potBound->setObjectName(QStringLiteral("lbl_potBound"));

        horizontalLayout_6->addWidget(lbl_potBound);

        lineEdit_potbound = new QLineEdit(centralWidget);
        lineEdit_potbound->setObjectName(QStringLiteral("lineEdit_potbound"));
        sizePolicy.setHeightForWidth(lineEdit_potbound->sizePolicy().hasHeightForWidth());
        lineEdit_potbound->setSizePolicy(sizePolicy);
        lineEdit_potbound->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_6->addWidget(lineEdit_potbound);


        gridLayout->addLayout(horizontalLayout_6, 1, 0, 1, 1);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(6);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        lbl_E0 = new QLabel(centralWidget);
        lbl_E0->setObjectName(QStringLiteral("lbl_E0"));

        horizontalLayout_10->addWidget(lbl_E0);

        lineEdit_E0 = new QLineEdit(centralWidget);
        lineEdit_E0->setObjectName(QStringLiteral("lineEdit_E0"));
        sizePolicy.setHeightForWidth(lineEdit_E0->sizePolicy().hasHeightForWidth());
        lineEdit_E0->setSizePolicy(sizePolicy);
        lineEdit_E0->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_10->addWidget(lineEdit_E0);


        gridLayout->addLayout(horizontalLayout_10, 1, 1, 1, 1);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setSpacing(6);
        horizontalLayout_12->setObjectName(QStringLiteral("horizontalLayout_12"));
        label_alphaBeamMax = new QLabel(centralWidget);
        label_alphaBeamMax->setObjectName(QStringLiteral("label_alphaBeamMax"));

        horizontalLayout_12->addWidget(label_alphaBeamMax);

        lineEdit_alphaBeamMax = new QLineEdit(centralWidget);
        lineEdit_alphaBeamMax->setObjectName(QStringLiteral("lineEdit_alphaBeamMax"));
        sizePolicy.setHeightForWidth(lineEdit_alphaBeamMax->sizePolicy().hasHeightForWidth());
        lineEdit_alphaBeamMax->setSizePolicy(sizePolicy);
        lineEdit_alphaBeamMax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_12->addWidget(lineEdit_alphaBeamMax);


        gridLayout->addLayout(horizontalLayout_12, 2, 0, 1, 1);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        lbl_numgpus = new QLabel(centralWidget);
        lbl_numgpus->setObjectName(QStringLiteral("lbl_numgpus"));

        horizontalLayout_9->addWidget(lbl_numgpus);

        spinBox_numGPUs = new QSpinBox(centralWidget);
        spinBox_numGPUs->setObjectName(QStringLiteral("spinBox_numGPUs"));
        spinBox_numGPUs->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_9->addWidget(spinBox_numGPUs);


        gridLayout->addLayout(horizontalLayout_9, 2, 1, 1, 1);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setSpacing(6);
        horizontalLayout_13->setObjectName(QStringLiteral("horizontalLayout_13"));
        label_sliceThickness = new QLabel(centralWidget);
        label_sliceThickness->setObjectName(QStringLiteral("label_sliceThickness"));

        horizontalLayout_13->addWidget(label_sliceThickness);

        lineEdit_sliceThickness = new QLineEdit(centralWidget);
        lineEdit_sliceThickness->setObjectName(QStringLiteral("lineEdit_sliceThickness"));
        sizePolicy.setHeightForWidth(lineEdit_sliceThickness->sizePolicy().hasHeightForWidth());
        lineEdit_sliceThickness->setSizePolicy(sizePolicy);
        lineEdit_sliceThickness->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_13->addWidget(lineEdit_sliceThickness);


        gridLayout->addLayout(horizontalLayout_13, 3, 0, 1, 1);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        lbl_numthreads = new QLabel(centralWidget);
        lbl_numthreads->setObjectName(QStringLiteral("lbl_numthreads"));

        horizontalLayout_8->addWidget(lbl_numthreads);

        spinBox_numThreads = new QSpinBox(centralWidget);
        spinBox_numThreads->setObjectName(QStringLiteral("spinBox_numThreads"));
        spinBox_numThreads->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_8->addWidget(spinBox_numThreads);


        gridLayout->addLayout(horizontalLayout_8, 3, 1, 1, 1);


        verticalLayout_2->addLayout(gridLayout);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setSpacing(6);
        horizontalLayout_11->setObjectName(QStringLiteral("horizontalLayout_11"));
        label_cellDim = new QLabel(centralWidget);
        label_cellDim->setObjectName(QStringLiteral("label_cellDim"));

        horizontalLayout_11->addWidget(label_cellDim);

        lineEdit_cellDimX = new QLineEdit(centralWidget);
        lineEdit_cellDimX->setObjectName(QStringLiteral("lineEdit_cellDimX"));
        sizePolicy.setHeightForWidth(lineEdit_cellDimX->sizePolicy().hasHeightForWidth());
        lineEdit_cellDimX->setSizePolicy(sizePolicy);
        lineEdit_cellDimX->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_11->addWidget(lineEdit_cellDimX);

        lineEdit_cellDimY = new QLineEdit(centralWidget);
        lineEdit_cellDimY->setObjectName(QStringLiteral("lineEdit_cellDimY"));
        sizePolicy.setHeightForWidth(lineEdit_cellDimY->sizePolicy().hasHeightForWidth());
        lineEdit_cellDimY->setSizePolicy(sizePolicy);
        lineEdit_cellDimY->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_11->addWidget(lineEdit_cellDimY);

        lineEdit_cellDimZ = new QLineEdit(centralWidget);
        lineEdit_cellDimZ->setObjectName(QStringLiteral("lineEdit_cellDimZ"));
        sizePolicy.setHeightForWidth(lineEdit_cellDimZ->sizePolicy().hasHeightForWidth());
        lineEdit_cellDimZ->setSizePolicy(sizePolicy);
        lineEdit_cellDimZ->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_11->addWidget(lineEdit_cellDimZ);


        verticalLayout_2->addLayout(horizontalLayout_11);


        horizontalLayout_4->addLayout(verticalLayout_2);

        btn_go = new QPushButton(centralWidget);
        btn_go->setObjectName(QStringLiteral("btn_go"));
        btn_go->setMinimumSize(QSize(0, 145));
        QFont font;
        font.setPointSize(30);
        font.setBold(true);
        font.setWeight(75);
        btn_go->setFont(font);

        horizontalLayout_4->addWidget(btn_go);


        verticalLayout->addLayout(horizontalLayout_4);


        verticalLayout_3->addLayout(verticalLayout);

        PRISMMainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(PRISMMainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 778, 22));
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
        PRISMMainWindow->setWindowTitle(QApplication::translate("PRISMMainWindow", "PRISMMainWindow", 0));
        lbl_atomsfile->setText(QApplication::translate("PRISMMainWindow", "Atoms File ", 0));
        btn_atomsfile_browse->setText(QApplication::translate("PRISMMainWindow", "Browse", 0));
        lbl_outputfile->setText(QApplication::translate("PRISMMainWindow", "Output File", 0));
        btn_outputfile_browse->setText(QApplication::translate("PRISMMainWindow", "Browse", 0));
        label->setText(QApplication::translate("PRISMMainWindow", "Interpolation factor (f)", 0));
        lbl_pixelSize->setText(QApplication::translate("PRISMMainWindow", "Pixel Size", 0));
        lineedit_pixelSize->setText(QApplication::translate("PRISMMainWindow", "0.1", 0));
        lbl_numfp->setText(QApplication::translate("PRISMMainWindow", "Number of Frozen Phonons", 0));
        lbl_potBound->setText(QApplication::translate("PRISMMainWindow", "Potential Bound", 0));
        lineEdit_potbound->setText(QApplication::translate("PRISMMainWindow", "1", 0));
        lbl_E0->setText(QApplication::translate("PRISMMainWindow", "Energy (keV)", 0));
        lineEdit_E0->setText(QApplication::translate("PRISMMainWindow", "80", 0));
        label_alphaBeamMax->setText(QApplication::translate("PRISMMainWindow", "Max Angle (mrad)", 0));
        lineEdit_alphaBeamMax->setText(QApplication::translate("PRISMMainWindow", "24", 0));
        lbl_numgpus->setText(QApplication::translate("PRISMMainWindow", "Num GPUs", 0));
        label_sliceThickness->setText(QApplication::translate("PRISMMainWindow", "Slice Thickness", 0));
        lineEdit_sliceThickness->setText(QApplication::translate("PRISMMainWindow", "2", 0));
        lbl_numthreads->setText(QApplication::translate("PRISMMainWindow", "Num Threads", 0));
        label_cellDim->setText(QApplication::translate("PRISMMainWindow", "Cell Dimension", 0));
        lineEdit_cellDimX->setText(QApplication::translate("PRISMMainWindow", "100", 0));
        lineEdit_cellDimY->setText(QApplication::translate("PRISMMainWindow", "100", 0));
        lineEdit_cellDimZ->setText(QApplication::translate("PRISMMainWindow", "80", 0));
        btn_go->setText(QApplication::translate("PRISMMainWindow", "Calculate!", 0));
    } // retranslateUi

};

namespace Ui {
    class PRISMMainWindow: public Ui_PRISMMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PRISMMAINWINDOW_H
