#!/bin/bash

# 1st param: MPI -n and C program first parameter
# 2nd param: MPI --hostfile
# 3rd param: {file}.c
# Remaining params: C program remaining parameters
echo "mpicc -fopenmp $3.c -o $3.o -lm"
mpicc -fopenmp $3.c -o $3.o -lm
echo "mpirun -n $1 --hostfile $2 $3.o $1 ${@:4}"
mpirun -n $1 --hostfile $2 $3.o $1 ${@:4}