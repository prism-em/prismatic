/********************************************************************************
** Form generated from reading UI file 'saveatomiccoordinatesdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SAVEATOMICCOORDINATESDIALOG_H
#define UI_SAVEATOMICCOORDINATESDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_SaveAtomicCoordinatesDialog
{
public:
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout;
    QVBoxLayout *verticalLayout_2;
    QLabel *label_2;
    QLabel *label;
    QVBoxLayout *verticalLayout;
    QLineEdit *lineEdit_Filename;
    QLineEdit *lineEdit_Comment;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *SaveAtomicCoordinatesDialog)
    {
        if (SaveAtomicCoordinatesDialog->objectName().isEmpty())
            SaveAtomicCoordinatesDialog->setObjectName(QStringLiteral("SaveAtomicCoordinatesDialog"));
        SaveAtomicCoordinatesDialog->resize(660, 157);
        verticalLayout_3 = new QVBoxLayout(SaveAtomicCoordinatesDialog);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        label_2 = new QLabel(SaveAtomicCoordinatesDialog);
        label_2->setObjectName(QStringLiteral("label_2"));

        verticalLayout_2->addWidget(label_2);

        label = new QLabel(SaveAtomicCoordinatesDialog);
        label->setObjectName(QStringLiteral("label"));

        verticalLayout_2->addWidget(label);


        horizontalLayout->addLayout(verticalLayout_2);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        lineEdit_Filename = new QLineEdit(SaveAtomicCoordinatesDialog);
        lineEdit_Filename->setObjectName(QStringLiteral("lineEdit_Filename"));

        verticalLayout->addWidget(lineEdit_Filename);

        lineEdit_Comment = new QLineEdit(SaveAtomicCoordinatesDialog);
        lineEdit_Comment->setObjectName(QStringLiteral("lineEdit_Comment"));

        verticalLayout->addWidget(lineEdit_Comment);


        horizontalLayout->addLayout(verticalLayout);


        verticalLayout_3->addLayout(horizontalLayout);

        buttonBox = new QDialogButtonBox(SaveAtomicCoordinatesDialog);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout_3->addWidget(buttonBox);


        retranslateUi(SaveAtomicCoordinatesDialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), SaveAtomicCoordinatesDialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), SaveAtomicCoordinatesDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(SaveAtomicCoordinatesDialog);
    } // setupUi

    void retranslateUi(QDialog *SaveAtomicCoordinatesDialog)
    {
        SaveAtomicCoordinatesDialog->setWindowTitle(QApplication::translate("SaveAtomicCoordinatesDialog", "Save Tiled Atomic Coordinates", Q_NULLPTR));
        label_2->setText(QApplication::translate("SaveAtomicCoordinatesDialog", "Filename", Q_NULLPTR));
        label->setText(QApplication::translate("SaveAtomicCoordinatesDialog", "Comment", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class SaveAtomicCoordinatesDialog: public Ui_SaveAtomicCoordinatesDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SAVEATOMICCOORDINATESDIALOG_H
