#!/bin/bash

set -e

echo "====== Building AMRaCut-Testbench ======"

mkdir -p build

cmake -DENABLE_VTK_FEATURES=ON -S ./testbench -B build
cmake --build build -- VERBOSE=1

echo "====== AMRaCut-Testbench build done ======"