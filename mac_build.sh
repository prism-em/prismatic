#! /bin/bash
export QT_PATH='/Users/ajpryor/Qt/5.9/'
mkdir prismatic-gui.app/Contents/Frameworks
mkdir prismatic-gui.app/Contents/PlugIns
mkdir prismatic-gui.app/Contents/PlugIns/platforms

cp -R ${QT_PATH}clang_64/lib/QtCore.framework prismatic-gui.app/Contents/Frameworks/
cp -R ${QT_PATH}clang_64/lib/QtGui.framework prismatic-gui.app/Contents/Frameworks/
cp -R ${QT_PATH}clang_64/lib/QtWidgets.framework prismatic-gui.app/Contents/Frameworks/
cp -R ${QT_PATH}clang_64/lib/QtPrintSupport.framework prismatic-gui.app/Contents/Frameworks/
# cp -R ${QT_PATH}/clang_64/plugins/platforms/libqcocoa.dylib prismatic-gui.app/Contents/PlugIns
cp -R ${QT_PATH}/clang_64/plugins/platforms/libqcocoa.dylib prismatic-gui.app/Contents/PlugIns/platforms/

echo "@[Paths]
Plugins = PlugIns" > prismatic-gui.app/Contents/Resources/qt.conf

install_name_tool -change @rpath/QtGui.framework/Versions/5/QtGui \
 @executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui \
  prismatic-gui.app/Contents/PlugIns/platforms/libqcocoa.dylib

install_name_tool -change @rpath/QtCore.framework/Versions/5/QtCore \
 @executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore \
  prismatic-gui.app/Contents/PlugIns/platforms/libqcocoa.dylib

install_name_tool -change @rpath/QtWidgets.framework/Versions/5/QtWidgets \
 @executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets \
  prismatic-gui.app/Contents/PlugIns/platforms/libqcocoa.dylib

install_name_tool -change @rpath/QtPrintSupport.framework/Versions/5/QtPrintSupport \
 @executable_path/../Frameworks/QtPrintSupport.framework/Versions/5/QtPrintSupport \
  prismatic-gui.app/Contents/PlugIns/platforms/libqcocoa.dylib

install_name_tool -id @executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore prismatic-gui.app/Contents/Frameworks/QtCore.framework/Versions/5/QtCore
install_name_tool -id @executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui prismatic-gui.app/Contents/Frameworks/QtGui.framework/Versions/5/QtGui
install_name_tool -id @executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets prismatic-gui.app/Contents/Frameworks/QtWidgets.framework/Versions/5/QtWidgets

install_name_tool -change @rpath/QtCore.framework/Versions/5/QtCore \
 @executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore \
  prismatic-gui.app/Contents/MacOS/prismatic-gui

install_name_tool -change @rpath/QtGui.framework/Versions/5/QtGui \
 @executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui \
  prismatic-gui.app/Contents/MacOS/prismatic-gui

install_name_tool -change @rpath/QtWidgets.framework/Versions/5/QtWidgets \
 @executable_path/../Frameworks/QtWidgets.framework/Versions/5/QtWidgets \
  prismatic-gui.app/Contents/MacOS/prismatic-gui

install_name_tool -change @rpath/QtCore.framework/Versions/5/QtCore \
  @executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore \
  prismatic-gui.app/Contents/Frameworks/QtGui.framework/Versions/5/QtGui

install_name_tool -change @rpath/QtGui.framework/Versions/5/QtGui \
  @executable_path/../Frameworks/QtGui.framework/Versions/5/QtGui \
  prismatic-gui.app/Contents/Frameworks/QtWidgets.framework/Versions/5/QtWidgets

install_name_tool -change @rpath/QtCore.framework/Versions/5/QtCore \
  @executable_path/../Frameworks/QtCore.framework/Versions/5/QtCore \
  prismatic-gui.app/Contents/Frameworks/QtWidgets.framework/Versions/5/QtWidgets