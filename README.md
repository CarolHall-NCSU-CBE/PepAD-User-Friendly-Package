# PepAD-User-Friendly-Package
PepAD: a user-friendly Monte Carlo algorithm for designing amyloid-forming self-assembling peptides
## Introduction
The Peptide Assembly Design (PepAD) is a Monte Carlo (MC)-based algorithm for discovering self-assembling peptides [1-3]. PepAD can design peptides that can self-assemble into the user-defined reference structures. Users can customize designs according to their experimental goals, such as excluding residues from mutations, restricting the number of certain amino acids, and designating NMR-constrained sites.

## Workflow
  A simplified flow chart describing the PepAD algorithm is shown in Fig 1. The algorithm begins with an initial structure peptide fibril structure. This structure typically consists of two stacked β-sheets and can be classified as a certain class of cross β-spine, defined by Sawaya and Eisenberg [4-5]. To initiate the design process, a random sequence is draped upon the backbone (-NH-Cα-CO-) of the initial structure. The binding affinity of the resulting fibril structure is evaluated using a score function (Γ<sub>score</sub>). Then one of the three types of MC moves are applied to search for new self-assembling peptides. After each move, the Γ<sub>score</sub> is evaluated for the new configuration (after sequence change move or sheet position perturbation move). The Monte-Carlo Metropolis criterion is then applied to accept or reject this new configuration. The entire process is iterated for multiple steps and the peptides with lowest Γ<sub>score</sub> are usually considered as the best design throughout the entire run.
<p align="center">
  <img width="300" alt="Flow chart of PepAD algorithm" src="https://github.com/user-attachments/assets/572abb95-0b87-435f-b154-17b3910b86e4">
  </p>
<p align="center"><b>Fig. 1.</b> Simplified flow chart describing the PepAD algorithm.</figcaption></p>

### Monte-Carlo Move  
  Three types of MC moves are applied to search for new self-assembling peptides: residue mutation, residue exchange, and sheet position perturbation. (1) Residue mutation: a random amino acid in the sequence (for all peptides) is mutated to another amino acid of same hydration type, (e.g. ALA is mutated to VAL (hydrophobic), SER is mutated to ASN (polar), or ARG is mutated to LYS (charged)); (2) Residue exchange: two random residues in the sequence are exchanged. (3) Sheet position perturbation move: the position of one β-sheet shifts in a random direction with a random displacement, while the position of another β-sheet is fixed.
<p align="center"> 
  <img width="500" alt="image" src="https://github.com/user-attachments/assets/e44bc028-d5a0-4896-b5cc-6f442344231e" />
  </p>
<p align="center"><b>Fig. 2.</b>  MC moves in PepAD.</figcaption></p>

### Score Function
In PepAD, after each MC move, a score function Γ<sub>score</sub> is used to quantify the binding affinity of the new sequence Γ<sub>score</sub> contains (i) a binding free energy term ΔG<sub>bind</sub> and (ii) an aggregation propensity term P<sub>agg</sub>. A weighting factor $$\lambda$$ is used to balance the contribution of ΔG<sub>bind</sub> and P<sub>agg</sub> to the Γ<sub>score</sub>.

$$ Γ_{score} = ΔG_{bind} - \lambda \times P_{agg} $$

The highly negative value of Γ<sub>score</sub> represents that the evaluated peptide sequence has a strong binding affinity at a given configuration. ΔG<sub>bind</sub> is calculated using the molecular mechanics-generalized Born surface area (MMGBSA) method [6–8]. P<sub>agg</sub> accounts for the intrinsic aggregation propensity of a sequence based on the Zyggregator method proposed by Dobson and Vendruscolo [9–11].

## Usage
### Compilation
The coding language is Fortran90. Users need to load Intel Fortran Compiler (ifx) for compilation. In the command line interface, navigate to the `/src/`, and use the following command:
 
>
    module load PrgEnv-intel         # load Intel compiler

    ifx -o PepAD main.f90            #  compile PepAD code


> **Note:** Code must be compiled in the same directory as `/lib/`, since parameter files are required at runtime. After compilation, the PepAD can be run from any directory.

### Prepare Input Files for PepAD
PepAD requires an `input.txt` file and an `initial_structure.pdb` to start searching. 

#### initial structure
Users must provide an initial structure PDB file for PepAD. The PDB file format must match that in `/example/`. The structure usually consists of inidividual peptides in a configuration of cross-β spine. Two recommended approaches for preparing the initial structure are as follows: 

(a) Search for existing amyloid fibril structures in the Protein Data Bank ([RCSB PDB: Homepage](https://www.rcsb.org/)). Users need to modify the format of the PDB file to match PepAD's required input format (see Fig. 3, left).

(b) Build artificial amyloid backbones using the peptide-building tool provided in this package, or prepare with other molecular modelling tools (UCSF Chimera or Packmol) (see Fig. 3, right).

<p align="center"> 
  <img width="420" alt="image" src="https://github.com/user-attachments/assets/3059ea8b-0e86-44dd-858b-50d21276b000" />
  <img width="420" alt="image" src="https://github.com/user-attachments/assets/e5b3ac63-57e5-465e-ac81-ac454022c324" />
  </p>
<p align="center"><b>Fig. 3. (Left) A 7-mer fibril structure built using crystal structure of peptide GNNQQNY (PDB ID: 2omm) and molecular dynamics simulation. (Right) A 14-mer fibril strucutre built using peptide-building tool in PepAD package.</b>  .</figcaption></p>

#### input.txt
The 'input.txt' contains the necessary parameters for designing a fibril structure.
**Table 1:** Input paremeters for PepAD
| Item | Description |
|------|-------------|
| AA per chain | Total number of amino acids in a single peptide chain |
| total chains | Total number of peptide chains |
| total sheets | Total number of β-sheets |
| initial_structure | Name of the initial structure PDB file |
| recalculation switch | Switch for restarting PepAD runs:<br>=0 → fresh run from step_start to step_end<br>=1 → restart from where the previous run ends |
| step start | The starting step number |
| step end | The ending step number. Thus, total steps = (step_end - step_start) + 1 |
| seed switch | Switch for employing a random seed:<br>=0 → use random seed based on local time<br>=1 → use provided random seed |
| random seed | The user-provided random seed when seed_switch = 1 |
| ekt_seq | Temperature factor $$k_{B}T$$ used in the Metropolis criterion for residue mutation/exchange |
| sheet switch | Probability of employing a sheet position perturbation move (0–1).<br>=0 → no sheet move<br>=0.6 → 60% chance of sheet move |
| ekt_sheet | Temperature factor $$k_{B}T$$ used in the Metropolis criterion for sheet position perturbation move |
| sheetmove interval | Number of steps to only perform residue mutation/exchange moves after a successful sheet position perturbation move |
| rmsd_max_x, rmsd_max_y, rmsd_max_z | Maximum displacement (absolute coordinates) a β-sheet can shift in x-, y-, or z-direction |
| dx, dy, dz | Maximum single-step displacement a β-sheet can shift in x-, y-, or z-direction |
| propensity weighting factor | λ, the weighting factor for the contribution of aggregation propensity in the score function |
| N_hydrophobic, N_polar, N_charged, N_other | Number of hydrophobic, polar, charged, and other residues in the desired sequence (excluding NMR-constrained residues) |
| Number of peptides in group (i) | Number of peptide chains in group (i) for energy calculation |
| chain ID of peptides in sheet (i) | Index numbers of peptide chains in group (i) |
| NMR constrained sites | Residues on these sites are constrained to the NMR amino acid pool. Example:<br>4, 5-6, 12 → sites 4, 5, 6, and 12 are constrained |
| Sites with no sequence change move | Residues on these sites are fixed to their initial amino acid throughout the run. Example:<br>1, 16 → amino acids on sites 1 and 16 remain unchanged |
| Must have amino acids | Sequence must include at least one of each listed amino acid (3-letter code, e.g., ALA, VAL). Enter `None` or leave blank for no restriction |
| NMR constrained amino acid pool | Define the set of amino acids for the NMR constrained sites (3-letter codes). Enter `None` or leave blank for no constraints |
| Restrictions on amino acids | Restrict the count of specific amino acids. Use 3-letter code + number. Example: `ASN3` (limit to 3 ASN). Enter `None` or leave blank for no restriction |

### Submission Script
Once the input files are prepared. User can submit the following example script to run PepAD:

>
	#!/bin/bash

    mkdir -p pdbfiles
    
    /path_to_PepAD/src/PepAD

### Output files
**Table 2:** PepAD output files 
|File extension| Description|
|----------------|------------|
|`Energy profile`| Records Steps, Sequences, Γ<sub>score</sub>, ΔG<sub>bind</sub>, ΔE<sub>bind</sub>, ΔTS<sub>conf</sub>, ΔTS<sub>conf</sub>, and P<sub>agg</sub>-I<sub>hydr</sub>|
|`Energy details`| Records Steps, Trials, Sequences, Γ<sub>score</sub>, ΔG<sub>bind</sub>, ΔE<sub>bind</sub>, ΔTS<sub>conf</sub>, ΔTS<sub>conf</sub>, and P<sub>agg</sub>-I<sub>hydr</sub>, MC moves, Trial results|
|`Minimal energy`| Records Steps, Trials, Γ<sub>score</sub>, and Sequences with minimal Γ<sub>score</sub> along the search|
|`RMSD profile`  | Records Steps, RMSD<sub>x</sub>, RMSD<sub>y</sub>, RMSD<sub>z</sub>, RMSD|
|`PDB files`     | PDB files for peptides with minimal scores during the evolution|
#### terms
- ΔG<sub>bind</sub>: Binding free energy

- ΔE<sub>bind</sub>: Binding energy

- ΔTS<sub>conf</sub>: Conformational entropy

- I<sub>hydr</sub>: Hydrophobic contribution in P<sub>agg</sub>

- P<sub>agg</sub>-I<sub>hydr</sub>: Aggregation propensity excepts the hydrophobic contribution

- RMSD<sub>x, y, z</sub>: Root Mean Square Displacement (RMSD) of a β-sheet in x, y, or z direction.

- RMSD: Total Root Mean Square Displacement (RMSD) of a β-sheet.



## Reference:
[1] S. Sarma, T.R. Sudarshan, V. Nguyen, A.S. Robang, X. Xiao, J.V. Le, M.E. Helmicki, A.K. Paravastu, C.K. Hall, Design of parallel 𝛽-sheet nanofibrils using Monte Carlo search, coarse-grained simulations, and experimental testing, Protein Science 33 (2024) e5102. https://doi.org/10.1002/pro.5102.

[2] X. Xiao, A.S. Robang, S. Sarma, J.V. Le, M.E. Helmicki, M.J. Lambert, R. Guerrero-Ferreira, J. Arboleda-Echavarria, A.K. Paravastu, C.K. Hall, Sequence patterns and signatures: Computational and experimental discovery of amyloid-forming peptides, PNAS Nexus 1 (2022) pgac263. https://doi.org/10.1093/pnasnexus/pgac263.

[3] X. Xiao, Y. Wang, D.T. Seroski, K.M. Wong, R. Liu, A.K. Paravastu, G.A. Hudalla, C.K. Hall, De novo design of peptides that coassemble into β sheet–based nanofibrils, Science Advances 7 (2021) eabf7668. https://doi.org/10.1126/sciadv.abf7668.

[4] M.R. Sawaya, S. Sambashivan, R. Nelson, M.I. Ivanova, S.A. Sievers, M.I. Apostol, M.J. Thompson, M. Balbirnie, J.J.W. Wiltzius, H.T. McFarlane, A.Ø. Madsen, C. Riekel, D. Eisenberg, Atomic structures of amyloid cross-β spines reveal varied steric zippers, Nature 447 (2007) 453–457. https://doi.org/10.1038/nature05695.

[5] R. Nelson, M.R. Sawaya, M. Balbirnie, A.Ø. Madsen, C. Riekel, R. Grothe, D. Eisenberg, Structure of the cross-β spine of amyloid-like fibrils, Nature 435 (2005) 773–778. https://doi.org/10.1038/nature03680.

[6] H. Gohlke, C. Kiel, D.A. Case, Insights into Protein–Protein Binding by Binding Free Energy Calculation and Free Energy Decomposition for the Ras–Raf and Ras–RalGDS Complexes, J. Mol. Biol. 330 (2003) 891–913. https://doi.org/10.1016/S0022-2836(03)00610-7.

[7]	G. Rastelli, A.D. Rio, G. Degliesposti, M. Sgobba, Fast and accurate predictions of binding free energies using MM-PBSA and MM-GBSA, J. Comput. Chem. 31 (2010) 797–810. https://doi.org/10.1002/jcc.21372.

[8]	X. Xiao, B. Zhao, P.F. Agris, C.K. Hall, Simulation study of the ability of a computationally‐designed peptide to recognize target tRNA Lys3 and other decoy tRNAs, Protein Sci. 25 (2016) 2243–2255. https://doi.org/10.1002/pro.3056.

[9] G.G. Tartaglia, A.P. Pawar, S. Campioni, C.M. Dobson, F. Chiti, M. Vendruscolo, Prediction of Aggregation-Prone Regions in Structured Proteins, J. Mol. Biol. 380 (2008) 425–436. https://doi.org/10.1016/j.jmb.2008.05.013.

[10]	G.G. Tartaglia, M. Vendruscolo, The Zyggregator method for predicting protein aggregation propensities, Chem. Soc. Rev. 37 (2008) 1395–1401. https://doi.org/10.1039/B706784B.

[11]	A.P. Pawar, K.F. DuBay, J. Zurdo, F. Chiti, M. Vendruscolo, C.M. Dobson, Prediction of “Aggregation-prone” and “Aggregation-susceptible” Regions in Proteins Associated with Neurodegenerative Diseases, J. Mol. Biol. 350 (2005) 379–392. https://doi.org/10.1016/j.jmb.2005.04.016.


