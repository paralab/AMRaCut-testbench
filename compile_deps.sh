#!/bin/bash

set -e

deps_install_dir="$PWD/install"

echo "====== Building and Installing AMRaCut ======"

cd AMRaCut

mkdir -p build

cmake -DAMRACUT_INTEGER_WIDTH=32 -DCMAKE_INSTALL_PREFIX="$deps_install_dir/AMRaCut" -S . -B build
cmake --build build --target install -- VERBOSE=1

echo "====== AMRaCut installation done ======"

echo "====== Building and Installing gmsh ======"

cd ../gmsh

mkdir -p build

cmake -DENABLE_FLTK=0 -DENABLE_BUILD_DYNAMIC=1 -DCMAKE_INSTALL_PREFIX="$deps_install_dir/gmsh" -S . -B build
cmake --build build --target install -- -j$(nproc) VERBOSE=1

echo "====== gmsh installation done ======"