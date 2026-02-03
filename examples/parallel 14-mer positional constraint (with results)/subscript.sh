#!/bin/bash
#BSUB -n 1
#BSUB -W 180:00
#BSUB -J PepAD
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -q hall

mkdir -p pdbfiles
./../src/PepAD



