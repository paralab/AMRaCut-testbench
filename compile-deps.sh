#!/bin/bash

set -e

deps_install_dir="$PWD/install"

cd AMRaCut

mkdir -p build

cmake -DAMRACUT_INTEGER_WIDTH=32 -DCMAKE_INSTALL_PREFIX="$deps_install_dir/AMRaCut" -S . -B build
cmake --build build --target install --verbose

cd ../gmsh
echo "TODO: gmsh compile"