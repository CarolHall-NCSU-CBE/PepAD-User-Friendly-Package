#!/usr/bin/env bash

set -e

source ~/.bashrc
conda activate pepad-tools # or conda activate path/to/pepad-tools

builder --seq AAAAAAAA --class 1 --dx 4.8 --dz 11.5 --x 1 \
    --chains 4 --cap 1 --core e --format 0 --output comp3.pdb
