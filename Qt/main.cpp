#include "prismmainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    PRISMMainWindow w;
    w.show();

    return a.exec();
}
