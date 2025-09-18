#!/bin/bash
#BSUB -n 1
#BSUB -W 10
#BSUB -q debug
#BSUB -o stdout.%J
#BSUB -e stderr.%J

module load intel/2017.1.132 intel_mpi/2017 PrgEnv-intel/2017.1.132
ifort -o PepAD main_v1.35-3.f90

