# jacobi-mpi

## Description

OpenMPI (v1.5) implementation of a simulation of thermal transmission in a 2D space.

OpenMP is also used to provide multithreading.

A singlethreaded and multithreaded implementation in Go can are available in [this repo](https://github.com/mcanalesmayo/jacobi-go).

Run the code by using the script ``comp_run.sh`` with the following parameters (ordered):
1. Number of MPI nodes
2. MPI hostfile
3. File containing the C code
4. Number of OMP threads for every node

## References

``jacobi_seq.c`` by **Rubén Gran Tejero** -- rgran@unizar.es

If you use this code for academic work, please also reference:
* **Rubén Gran Tejero** -- rgran@unizar.es
