#ifndef PRISMMAINWINDOW_H
#define PRISMMAINWINDOW_H

#include <QMainWindow>
#include <iostream>
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
    void printMessage(){
        std::cout << "Starting Calculation" << std::endl;
    }
private:
    Ui::PRISMMainWindow *ui;
};

#endif // PRISMMAINWINDOW_H
