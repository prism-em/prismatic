/********************************************************************************
** Form generated from reading UI file 'prism_progressbar.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PRISM_PROGRESSBAR_H
#define UI_PRISM_PROGRESSBAR_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_prism_progressbar
{
public:
    QHBoxLayout *horizontalLayout_3;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_2;
    QLabel *lbl_Algorithm;
    QLabel *lbl_Description;
    QProgressBar *progressBar;
    QHBoxLayout *horizontalLayout;
    QLabel *lbl_calcStatus;

    void setupUi(QDialog *prism_progressbar)
    {
        if (prism_progressbar->objectName().isEmpty())
            prism_progressbar->setObjectName(QStringLiteral("prism_progressbar"));
        prism_progressbar->resize(231, 122);
        horizontalLayout_3 = new QHBoxLayout(prism_progressbar);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        lbl_Algorithm = new QLabel(prism_progressbar);
        lbl_Algorithm->setObjectName(QStringLiteral("lbl_Algorithm"));

        horizontalLayout_2->addWidget(lbl_Algorithm);


        verticalLayout->addLayout(horizontalLayout_2);

        lbl_Description = new QLabel(prism_progressbar);
        lbl_Description->setObjectName(QStringLiteral("lbl_Description"));

        verticalLayout->addWidget(lbl_Description);

        progressBar = new QProgressBar(prism_progressbar);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        progressBar->setValue(24);

        verticalLayout->addWidget(progressBar);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        lbl_calcStatus = new QLabel(prism_progressbar);
        lbl_calcStatus->setObjectName(QStringLiteral("lbl_calcStatus"));

        horizontalLayout->addWidget(lbl_calcStatus);


        verticalLayout->addLayout(horizontalLayout);


        horizontalLayout_3->addLayout(verticalLayout);


        retranslateUi(prism_progressbar);

        QMetaObject::connectSlotsByName(prism_progressbar);
    } // setupUi

    void retranslateUi(QDialog *prism_progressbar)
    {
        prism_progressbar->setWindowTitle(QApplication::translate("prism_progressbar", "Dialog", Q_NULLPTR));
        lbl_Algorithm->setText(QApplication::translate("prism_progressbar", "PRISM:", Q_NULLPTR));
        lbl_Description->setText(QApplication::translate("prism_progressbar", "Computing Projected Potential", Q_NULLPTR));
        lbl_calcStatus->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class prism_progressbar: public Ui_prism_progressbar {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PRISM_PROGRESSBAR_H
