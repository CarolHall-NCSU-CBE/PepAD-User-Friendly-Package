#!/bin/bash
#BSUB -n 1
#BSUB -W 140:00
#BSUB -J PepADv1372
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -q hall

module load PrgEnv-intel
../../src/PepAD



