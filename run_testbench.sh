#!/bin/bash

set -e

export OMP_NUM_THREADS=1


mpirun -np 3 --oversubscribe ./build/testbench_main /home/budvin/research/Partitioning/mesh_generator/generated_tet_100x100x2.mesh
# mpirun -np 11 --oversubscribe ./build/testbench_main /home/budvin/research/Partitioning/Meshes/10k_hex/99898_sf_hexa.mesh
