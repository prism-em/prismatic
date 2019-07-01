#ifndef SAVEATOMICCOORDINATESDIALOG_H
#define SAVEATOMICCOORDINATESDIALOG_H

#include <QDialog>

namespace Ui
{
class SaveAtomicCoordinatesDialog;
}

class SaveAtomicCoordinatesDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SaveAtomicCoordinatesDialog(QWidget *parent = 0);
    ~SaveAtomicCoordinatesDialog();
public slots:
    void SaveAtomCoords();
    void setFilenameText(QString);
    void setCommentText(QString);
signals:
    void signalSaveAtomCoords(QString, QString);

private:
    Ui::SaveAtomicCoordinatesDialog *ui;
};

#endif // SAVEATOMICCOORDINATESDIALOG_H
