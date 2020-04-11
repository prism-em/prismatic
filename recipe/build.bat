setlocal EnableDelayedExpansion

if "%1"=="cli" (set cli=1 && set build=build_cli) else (set cli=0)

if "%1"=="gui" (set gui=1 && set build=build_gui) else (set gui=0)

echo cli is %cli% and gui is %gui%

set MENU_DIR=%PREFIX%\Menu
if not exist %MENU_DIR% mkdir %MENU_DIR%

if %gui%==1 (
  copy %RECIPE_DIR%\prismatic.ico %MENU_DIR%
  copy %RECIPE_DIR%\menu-windows.json %MENU_DIR%\prismatic_gui.json
)

:: Make a build folder and change to it.
mkdir %build% && cd %build%

:: Configure using the CMakeFiles
cmake -G "NMake Makefiles" ^
      -D PRISMATIC_ENABLE_CLI=%cli% ^
      -D PRISMATIC_ENABLE_GUI=%gui% ^
	  -D PRISMATIC_ENABLE_PYPRISMATIC=0 ^
      -D PRISMATIC_ENABLE_GPU=0 ^
      -D PRISMATIC_ENABLE_DOUBLE_PRECISION=0 ^
      -D HDF5_DIR=%LIBRARY_PREFIX%\cmake\hdf5 ^
      -D Qt5Widgets_DIR=%LIBRARY_PREFIX%\lib\cmake\Qt5Widgets ^
      -D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
      -D CMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
      -D CMAKE_BUILD_TYPE=Release ^
      %SRC_DIR%
if errorlevel 1 exit 1

:: Build!
nmake
if errorlevel 1 exit 1

:: Install!
nmake install
if errorlevel 1 exit 1

:: exit when don't build with double precision (only for cli)
if not %cli% == 1 exit 0

cd .. && mkdir build_double && cd build_double

cmake -G "NMake Makefiles" ^
      -D PRISMATIC_ENABLE_CLI=%cli% ^
      -D PRISMATIC_ENABLE_GUI=%gui% ^
      -D PRISMATIC_ENABLE_GPU=0 ^
      -D PRISMATIC_ENABLE_DOUBLE_PRECISION=1 ^
      -D HDF5_DIR=%LIBRARY_PREFIX%\cmake\hdf5 ^
      -D Qt5Widgets_DIR=%LIBRARY_PREFIX%\lib\cmake\Qt5Widgets ^
      -D OUTPUT_NAME="prismatic-double" ^
      -D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
      -D CMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
      -D CMAKE_BUILD_TYPE=Release ^
      %SRC_DIR%
if errorlevel 1 exit 1

:: Build!
nmake
if errorlevel 1 exit 1

:: Install!
nmake install
if errorlevel 1 exit 1
