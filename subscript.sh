#!/bin/bash
#BSUB -n 1
#BSUB -W 180:00
#BSUB -J PepAD
#BSUB -o stdout.%J
#BSUB -e stderr.%J


mkdir -p pdbfiles
./../src/PepAD_exe



