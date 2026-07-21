#!/bin/bash
# Example Intel Fortran compilation for the archived PepAD v1.37 source.
module load PrgEnv-intel
ifort -O2 -o PepAD main.f90
