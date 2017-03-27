#include <QApplication>
#include "prismmainwindow.h"
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    PRISMMainWindow w;
    w.setGeometry(0,0,950,650);
    w.show();
    return a.exec();
}
