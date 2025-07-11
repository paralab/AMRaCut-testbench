#!/bin/bash

set -e

VTK_INSTALL_DIR_PATH=/home/budvin/bin/VTK-9.3.0/build

echo "====== Building AMRaCut-Testbench ======"

mkdir -p build

cmake -DENABLE_VTK_FEATURES=ON -DVTK_INSTALL_DIR_PATH=${VTK_INSTALL_DIR_PATH} -S ./testbench -B build
cmake --build build

echo "====== AMRaCut-Testbench build done ======"