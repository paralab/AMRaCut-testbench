#!/bin/bash

set -e

echo "====== Building AMRaCut-Testbench ======"

mkdir -p build

cmake -S ./testbench -B build
cmake --build build -- VERBOSE=1

echo "====== AMRaCut-Testbench build done ======"