#include <QApplication>
#include "prismmainwindow.h"
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    PRISMMainWindow w;
   //
    w.show();
    w.setGeometry(0,0,950,650);
    return a.exec();
}
