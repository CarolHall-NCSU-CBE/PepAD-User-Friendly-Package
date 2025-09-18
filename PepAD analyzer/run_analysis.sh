#/bin/bash
conda activate pepad_analysis

# --top [n]: report "n" unique peptides with the lowest scores ranked in order, where "n" is an integer. 
# --plot [score, rmsd, both, none]: plot score versus step, RMSD versus step, or both versus step. Enter "none" for no plots.
# --score_rolling [i]: smooth the score plot with a rolling average over "i" steps.
# --rmsd_rolling [i]: smooth the rmsd plot with a rolling average over "i" steps.

python pepad_analyzer.py --top 10 --plot both --score_rolling 100 --rmsd_rolling 100