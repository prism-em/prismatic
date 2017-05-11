// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include <QApplication>
#include "prismmainwindow.h"
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    PRISMMainWindow w;
   //
    w.show();
//    w.setGeometry(100,100,950,650);
    w.setGeometry(100,100,950,750);

    return a.exec();
}
