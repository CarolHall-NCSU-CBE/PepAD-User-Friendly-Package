#!/usr/bin/env bash

set -e

source ~/.bashrc
conda activate pepad-tools # or conda activate path/to/pepad-tools

# --seq: peptide sequence in one-letter amino-acid codes
# --class: steric-zipper class from 1 to 10
# --dx: strand-strand distance along x in Å
# --dz: sheet-sheet distance along z in Å
# --x: sheet-2 shifts in x direction in unit of 0.5*dx
# --y: sheet-2 shifts in y direction in unit of residue spacing (3.465 Å)
# --chains: number of strands in each sheet
# --cap: 0 uncapped, 1 ACE+NME, or 2 ACE+NHE
# --core:  packing e = even number packed inside, o = odd number packed inside"
# --format: 0 PepAD PDB format; 1 AMBER format
# --output: output PDB filename


builder --seq GNNQQNY --class 1 --dx 4.8 --dz 11.0 --x -1 --y 1 \
    --chains 4 --cap 2 --core e --format 1 --output pep2.pdb
