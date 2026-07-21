#!/bin/bash
#BSUB -n 1
#BSUB -W 30
#BSUB -J env
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -q hall

module load conda
conda env create -f env.yml
