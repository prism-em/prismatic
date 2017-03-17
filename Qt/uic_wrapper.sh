#!/bin/sh
DYLD_LIBRARY_PATH=/Users/ajpryor/Qt5/qtbase/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH
QT_PLUGIN_PATH=/Users/ajpryor/Qt5/qtbase/plugins${QT_PLUGIN_PATH:+:$QT_PLUGIN_PATH}
export QT_PLUGIN_PATH
exec /Users/ajpryor/Qt5/qtbase/bin/uic "$@"
