
:: https://github.com/jkriege2/JKQtPlotter/issues/35
set "CXXFLAGS= -MD"
set "CMAKE_GENERATOR=NMake Makefiles"

%PYTHON% setup.py build_ext ^
		-DHDF5_DIR=%LIBRARY_PREFIX%\cmake\hdf5 ^
		-DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
		-DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
		install
if errorlevel 1 exit 1
