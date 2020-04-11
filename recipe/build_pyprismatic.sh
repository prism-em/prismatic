if [[ $1 == "gpu" ]]; then
	# remove -std=c++17 from CXXFLAGS for compatibility with nvcc
	echo "CXXFLAGS: $CXXFLAGS"
	export CXXFLAGS="$(echo $CXXFLAGS | sed -e 's/ -std=[^ ]*//')"
	echo "CXXFLAGS after removing -std=c++17: $CXXFLAGS"

	$PYTHON setup.py build_ext \
			-DCMAKE_PREFIX_PATH=$PREFIX \
			-DCMAKE_INSTALL_PREFIX=$PREFIX \
	  -DCUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
			install --enable-gpu
else
	$PYTHON setup.py build_ext \
			-DCMAKE_PREFIX_PATH=$PREFIX \
			-DCMAKE_INSTALL_PREFIX=$PREFIX \
			install
fi

