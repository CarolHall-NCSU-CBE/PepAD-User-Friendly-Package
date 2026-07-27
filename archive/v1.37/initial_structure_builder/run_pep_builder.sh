#!/bin/bash
#BSUB -n 1
#BSUB -W 10
#BSUB -J PepAD
#BSUB -o stdout.%J
#BSUB -e stderr.%J
#BSUB -q hall

source ~/.bashrc
conda activate /usr/local/usrapps/hall2/hwang226/env_buildpep

# -seq "peptide sequence"
# -c   "Fibril class"
# -sh  "the second sheet shifts in x-direction by -1, 0, or 1 residue"
# -n   "number of strands in a sheet"
# -p   "terminal patch (caps), 0 = no cap, 1 = ACE+NMR, 2 = ACE+NHE"
# -r   "residue packing e = even number packed inside, o = odd number packed inside"
# -f   "format flag 0 = PepAD, 1 = AMBER"
# -d1  "strand-strand distance (Angstrom)"
# -d2  "sheet-sheet distance (Angstrom)"
# -d3  "the second sheet shift in y-direction (along fibril axis) by a distance (Angstrom)"
# -o   "output file name"

python Initial_structure_builder.py -seq "GNNQQNY" \
    -c "1" -sh "-1" -n "4" -p "2" -r "e" -f "1" \
    -d1 "4.8" -d2 "10" -d3 "2.4" -o "pep2"

conda deactivate