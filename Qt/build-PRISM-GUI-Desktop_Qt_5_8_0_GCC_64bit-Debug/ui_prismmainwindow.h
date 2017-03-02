/********************************************************************************
** Form generated from reading UI file 'prismmainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PRISMMAINWINDOW_H
#define UI_PRISMMAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PRISMMainWindow
{
public:
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout_2;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *lbl_atomsfile;
    QLineEdit *lineedit_atomsfile;
    QPushButton *btn_atomsfile_browse;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QSpacerItem *horizontalSpacer;
    QLineEdit *lineedit_f;
    QHBoxLayout *horizontalLayout_3;
    QLabel *lbl_outputfile;
    QLineEdit *lineedit_outputfile;
    QPushButton *pushButton;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *PRISMMainWindow)
    {
        if (PRISMMainWindow->objectName().isEmpty())
            PRISMMainWindow->setObjectName(QStringLiteral("PRISMMainWindow"));
        PRISMMainWindow->resize(400, 300);
        centralWidget = new QWidget(PRISMMainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        verticalLayout_2 = new QVBoxLayout(centralWidget);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
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

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout_2->addWidget(label);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Maximum, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);

        lineedit_f = new QLineEdit(centralWidget);
        lineedit_f->setObjectName(QStringLiteral("lineedit_f"));

        horizontalLayout_2->addWidget(lineedit_f);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        lbl_outputfile = new QLabel(centralWidget);
        lbl_outputfile->setObjectName(QStringLiteral("lbl_outputfile"));

        horizontalLayout_3->addWidget(lbl_outputfile);

        lineedit_outputfile = new QLineEdit(centralWidget);
        lineedit_outputfile->setObjectName(QStringLiteral("lineedit_outputfile"));

        horizontalLayout_3->addWidget(lineedit_outputfile);


        verticalLayout->addLayout(horizontalLayout_3);

        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QStringLiteral("pushButton"));

        verticalLayout->addWidget(pushButton);


        verticalLayout_2->addLayout(verticalLayout);

        PRISMMainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(PRISMMainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 400, 22));
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
        lbl_atomsfile->setText(QApplication::translate("PRISMMainWindow", "Atoms File", Q_NULLPTR));
        btn_atomsfile_browse->setText(QApplication::translate("PRISMMainWindow", "Browse", Q_NULLPTR));
        label->setText(QApplication::translate("PRISMMainWindow", "iInterpolation factor (f)", Q_NULLPTR));
        lbl_outputfile->setText(QApplication::translate("PRISMMainWindow", "Output File", Q_NULLPTR));
        pushButton->setText(QApplication::translate("PRISMMainWindow", "PushButton", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class PRISMMainWindow: public Ui_PRISMMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PRISMMAINWINDOW_H
