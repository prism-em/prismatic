// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#define PRISM_BUILDING_GUI 1
#include <QApplication>
#include "prismmainwindow.h"
#include <QFontDatabase>
#include <QFont>
#include <QFile>
#include <QTextStream>

int main(int argc, char *argv[])
{
    QGuiApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QApplication a(argc, argv);
    QFontDatabase database;
    //not all fonts need to be loaded, probably
    int id1 = database.addApplicationFont(":/fonts/Roboto-Regular.ttf");
    int id2 = database.addApplicationFont(":/fonts/Roboto-Black.ttf");
    int id3 = database.addApplicationFont(":/fonts/Roboto-Light.ttf");
    int id4 = database.addApplicationFont(":/fonts/Roboto-Medium.ttf");
    QFont font = QFont("Roboto");
    a.setFont(font);
    PRISMMainWindow w;

    QFile file(":/dark.qss");
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream stream(&file);
    a.setStyleSheet(stream.readAll());

    //add a comment
    w.show();
    w.setGeometry(100, 100, 850, 700);

    return a.exec();
}
