#!/bin/bash

set -e

export OMP_NUM_THREADS=1


mpirun -np 12 --oversubscribe ./build/testbench_main /home/budvin/research/Partitioning/mesh_generator/generated_tet_100x100x2.mesh test.json
# mpirun -np 11 --oversubscribe ./build/testbench_main /home/budvin/research/Partitioning/Meshes/10k_hex/99898_sf_hexa.mesh test.json
