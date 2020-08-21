# setting the build flags
if [ "$1" = "" ] ; then
	echo "This script requires information on which mode to install."
	echo "usage: combined_build_install.sh [build opt]"
	echo "  [build opt] can be: gui or cli"
	exit 0 
fi
if [[ $1 == "cli" ]]; then
	enable_cli=1
	enable_gui=0
	build_path="build_cli"
elif [[ $1 == "gui" ]]; then
	enable_cli=0
	enable_gui=1
	build_path="build_gui"
else
	echo "option $1 is unrecognised"
	exit 0
fi 

if [[ $cuda_compiler_version != "None" && $OSTYPE == "linux-gnu" ]]; then
	enable_gpu=1
	# remove -std=c++17 from CXXFLAGS for compatibility with nvcc
	echo "CXXFLAGS: $CXXFLAGS"
	export CXXFLAGS="$(echo $CXXFLAGS | sed -e 's/ -std=[^ ]*//')"
	echo "CXXFLAGS after removing -std=c++17: $CXXFLAGS"
else
	enable_gpu=0
fi


# build process
mkdir $build_path && cd $build_path 

cmake -D PRISMATIC_ENABLE_GUI=$enable_gui \
	-D PRISMATIC_ENABLE_CLI=$enable_cli \
	-D PRISMATIC_ENABLE_GPU=$enable_gpu \
	-D PRISMATIC_ENABLE_PYPRISMATIC=0 \
	-D CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
	-D CMAKE_INSTALL_PREFIX=$PREFIX \
	-D CMAKE_PREFIX_PATH=${PREFIX} \
	../

make -j${CPU_COUNT}

make install

# exit when don't build with double precision (only for cli)
if [[ "$1" != "cli" ]] ; then exit 0 ; fi

echo "Building with double precision"

build_path="build_cli_double"
cd .. && mkdir $build_path  && cd $build_path 

cmake -D PRISMATIC_ENABLE_GUI=$enable_gui \
	-D PRISMATIC_ENABLE_CLI=$enable_cli \
	-D PRISMATIC_ENABLE_GPU=$enable_gpu \
	-D PRISMATIC_ENABLE_PYPRISMATIC=0 \
	-D CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
	-D PRISMATIC_ENABLE_DOUBLE_PRECISION=1 \
	-D OUTPUT_NAME="prismatic-double"\
	-D CMAKE_INSTALL_PREFIX=$PREFIX \
	-D CMAKE_PREFIX_PATH=${PREFIX} \
	../

make -j${CPU_COUNT}

make install
