#!/bin/bash

export CC=gcc
export CXX=g++

# specify the build directory
build_type=RELEASE
build_dir=$PWD/build
install_dir=${PWD}/program

# remove build_dir if it already exists
if [ -d $build_dir ] 
    then 
    rm -r $build_dir
fi

if [ -d $install_dir ]
    then 
    rm -rf $install_dir
fi

# make the build directory
mkdir -p $build_dir

# configure the build with cmake
cd $build_dir
cmake .. -DCMAKE_BUILD_TYPE=${build_type} \
         -DCMAKE_INSTALL_PREFIX=${install_dir}

# make with 8 threads
make -j 8 

make install