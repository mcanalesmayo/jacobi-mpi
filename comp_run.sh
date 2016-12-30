#!/bin/bash

# 1st param: MPI -n and first C programa parameter
# 2nd param: MPI --hostfile
# 3rd param: {file}.c
# Remaining params: C remaining program parameters
echo "mpicc $3.c -o $3.o -lm"
mpicc $3.c -o $3.o -lm
echo "mpirun -n $1 --hostfile $2 $3.o $1 ${@:4}"
mpirun -n $1 --hostfile $2 $3.o $1 ${@:4}