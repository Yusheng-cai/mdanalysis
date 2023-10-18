#!/bin/bash

export CC=gcc
export CXX=g++

# specify the build directory
build_type=RELEASE
build_dir=$PWD/${build_type}/
install_dir=${HOME}/programs/mdanalysis/
fftw3_dir=${HOME}/programs/fftw-3.3.10/

# remove build_dir if it already exists
[[ -d ${build_dir} ]] && rm -rf ${build_dir}
mkdir -p $build_dir
cd $build_dir

cmake .. \
	-DCMAKE_BUILD_TYPE=${build_type} \
	-DCMAKE_INSTALL_PREFIX=${install_dir} \
	-DFFTW3_DIR=${fftw3_dir}

# make with 8 threads
make -j 24
#make test
make install
