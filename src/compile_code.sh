#!/bin/bash
#BSUB -n 1
#BSUB -W 10
#BSUB -q debug
#BSUB -o stdout.%J
#BSUB -e stderr.%J

module load intel/2017.1.132 intel_mpi/2017 PrgEnv-intel/2017.1.132
ifort -O2 -o PepAD main_v1.42.f90
