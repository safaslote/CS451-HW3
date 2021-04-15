# CS451-HW3
Written by Safa Slote
This program will be using MPI on a Gauss Elimination algorthm

How to compile code:
mpicc -o gauss_mpi gauss_serial.c
mpirun -np <enter num of processors> gauss_mpi <random seed> <size of matrix>
