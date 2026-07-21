# PepAD v1.37 archive

> [!IMPORTANT]
> This directory preserves PepAD v1.37 for reproducibility. It uses the v1.37
> positional input format and is not the current release. For new calculations,
> return to the [repository README](../../README.md) and use PepAD v1.42.

# PepAD-User-Friendly-Package
## Introduction  
 
Peptide Assembly Design (PepAD) is a Monte Carlo (MC)–based algorithm for discovering amyloid-forming peptides that self-assemble into user-defined reference structures [1–3]. PepAD designs peptide by optimizing the binding affinity within a predefined fibril backbone, which is composed of two stacked β-sheets.

PepAD allows users to design peptides with specified sequence length, amino acid composition, and positional or compositional constraints. During the design process, residues on the peptide backbone are modified through MC moves. The new sequence is evaluated using a physics-based score function that quantifies the peptide self-assembly propensity.

To run PepAD, users first need to compile the code in HPC with Intel@ Fortran compiler. Then need to provide an initial fibril structure (PDB format) and a text-based input file that defines parameters and constraints. PepAD outputs peptide sequences with their corresponding scores, as well as structural files for the top-scoring designs, which can be further used in downstream simulations or experiments.


## PepAD algorithm
  PepAD starts from an initial peptide fibril structure, typically consisting of two stacked β-sheets corresponding to a specific class of cross-β spine[4-5]. At beign, a random sequence is draped upon the backbone (–NH–Cα–CO–), and its binding affinity is evaluated using a score function (Γ<sub>score</sub>). PepAD then explores new sequences using three MC moves—residue mutation, residue exchange, and sheet position perturbation. After each move, the Γ<sub>score</sub> is evaluated for the new configuration, and the new design is accepted or rejected based on Metropolis criterion. This process is iterated over multiple steps, and peptides with the lowest Γ<sub>score</sub> are identified as the best designs.

<p align="center">
  <img width="450" alt="image" src="https://github.com/user-attachments/assets/b0a786f8-a65d-48dd-a610-cdfaaca9e7fb" />
  </p>
<p align="center"><b>Fig. 1.</b> Simplified flow chart describing the PepAD algorithm.</figcaption></p>

### Monte-Carlo Move  
Three types of MC moves (Fig. 1 Perform one MC move) are applied to search for new self-assembling peptides: residue mutation, residue exchange, and sheet position perturbation. 

(1) Residue mutation: a random amino acid in the sequence (for all peptides) is mutated to another amino acid of same hydration type, (e.g. ALA is mutated to VAL (hydrophobic), SER is mutated to ASN (polar), or ARG is mutated to LYS (charged)); 

(2) Residue exchange: two random residues in the sequence are exchanged. 

(3) Sheet position perturbation move: the position of one β-sheet shifts in a random direction with a random displacement, while the position of another β-sheet is fixed.


### Score Function
In PepAD, after each MC move, a score function Γ<sub>score</sub> is used to quantify the binding affinity of the new sequence Γ<sub>score</sub> contains (i) a binding free energy term ΔG<sub>binding</sub> and (ii) an aggregation propensity term P<sub>aggregation</sub>. A weighting factor $$\lambda$$ is used to balance the contribution of ΔG<sub>binding</sub> and P<sub>aggregation</sub> to the Γ<sub>score</sub>.

$$ Γ_{score} = ΔG_{binding} - \lambda \times P_{aggregation} $$

The highly negative value of Γ<sub>score</sub> represents that the evaluated peptide sequence has a strong binding affinity at a given configuration. ΔG<sub>binding</sub> is calculated using the molecular mechanics-generalized Born surface area (MMGBSA) method [6–8]. P<sub>aggregation</sub> accounts for the intrinsic aggregation propensity of a sequence based on the Zyggregator method proposed by Dobson and Vendruscolo [9–11].

## Usage
<p align="center">
 <img width="1200"  alt="image" src="https://github.com/user-attachments/assets/982d5573-7dac-4381-a9ce-651992b16658" />
 </p>
<p align="center"><b>Fig. 2.</b> A workflow chart describing how to design peptides.</figcaption></p>


### Compile code
The coding language is Fortran90. Users need to load Intel Fortran Compiler (ifx) for compilation. In the command line interface, navigate to the `/src/`, and use the following command:
 
>

	# Example compilation using the Intel Fortran compiler on a Linux-based system
	
	# Users should modify these commands according to their local environment. 

	module load PrgEnv-intel         # load Intel compiler

    ifx -o PepAD main_v1.37.f90      # compile PepAD code


> **Note:** Code must be compiled in the same directory as `/lib/`, since parameter files are required at runtime. After compilation, the PepAD can be run from any directory.

### Set up input files
PepAD requires two files

(1) `initial structure PDB` file

(2) `input.txt` file 

#### initial structure
Users must provide an initial structure PDB file for PepAD. The PDB file format must match that in `/example/`. The structure usually consists of inidividual peptides in a configuration of cross-β spine. Two recommended approaches for preparing the initial structure are as follows: 

(a) Search for existing amyloid fibril structures in the Protein Data Bank ([RCSB PDB: Homepage](https://www.rcsb.org/)) (Fig 3). For example, the Sup35 prion segment GNNQQNY structure (PDB ID: 2OMM) is replicated to generate a larger Class-1 cross-β fibril model. The resulting structure was relaxed using explicit-solvent atomistic MD simulation to remove atomic overlaps. The relaxed structure can be further used as the initial structure for PepAD. Users also need to modify the format of the PDB file to match PepAD's required input format (see Fig. 3, top).

<p align="center">
  <img width="650" alt="image" src="https://github.com/user-attachments/assets/c476b3f8-90dd-411e-837f-1fc97d476770" />
  </p>
<p align="center"><b>Fig. 3.</b> Prepare initial structure using existing fibril structure from Protein Data Bank.</figcaption></p>


(b) Build artificial amyloid backbones using molecular modeling tools such as UCSF Chimera and Packmol (Fig. 4). For example, Aβ(16–22) KLVFFAE was built with UCSF Chimera, packed into a two-layer Class-8 antiparallel β-sheet fibril using Packmol, and relaxed by explicit-solvent MD simulation. To reduce PepAD computational cost, the middle eight peptides were extracted from the relaxed fibril and saved as comp1.pdb. The amyloid backbone can also be generated using the Initial Structure Builder provided in this package.

<p align="center"> 
  <img width="800" alt="image" src="https://github.com/user-attachments/assets/3b1b6026-a947-4a42-98e6-0e8600c831d6" />
  </p>
<p align="center"><b>Fig. 4.</b> Prepare initial structure using molecular modeling tool.</figcaption></p>


#### input.txt
The 'input.txt' contains the necessary parameters for designing a fibril structure.
**Table 1:** Input paremeters for PepAD
|  Panel |  Item | Description |
|--------|-------|-------------|
| Initial structure info | Number of residues per chain (including caps) | Total number of residues per peptide chain, including terminal caps if present. |
| Initial structure info | Number of chains in total | Total number of peptide chains in the initial structure PDB file. |
| Initial structure info | Input PDB file name | Name of the initial structure PDB file. |
| Initial structure info | Restart | Run mode selection: <br>`0` = fresh run<br> `1` = restart from the last completed step<br> |
| Initial structure info | Group 1 | Indices of peptide chains assigned to group 1 for binding free energy calculation. |
| Initial structure info | Group 2 | Indices of peptide chains assigned to group 2 for binding free energy calculation. |
| MC info | Number of MC moves | Integer. Total number of MC steps. |
| MC info | kBT(sequence) | Temperature factor `kBT` used in the Metropolis criterion for residue mutation and residue exchange moves. |
| MC info | kBT(sheet move) | Temperature factor `kBT` used in the Metropolis criterion for sheet position perturbation moves. |
| MC info | Allowed maximum RMSD (Å) | Maximum allowed displacement of one β-sheet relative to the other β-sheet during sheet position perturbation. |
| MC info | Aggregation propensity weight | Weighting factor `λ` that controls the contribution of aggregation propensity to the score function. |
| Compositions | Hydrophobic | Number of hydrophobic residues participating in sequence-change moves: `GLY`, `ALA`, `LEU`, `VAL`, `ILE`, `MET`, `PHE`, `TYR`, and `TRP`. |
| Compositions | Polar | Number of polar residues participating in sequence-change moves: `ASN`, `GLN`, `SER`, `THR`, and `HIE`. |
| Compositions | Charged | Number of charged residues participating in sequence-change moves: `ARG`, `LYS`, `GLU`, and `ASP`. |
| Compositions | Other | Number of residues participating in sequence-change moves that are not classified as hydrophobic, polar, or charged: `CYS` and `PRO`. |
| Constraints | Single-site positional constraint | Sequence sites where residues remain unchanged entire PepAD run. Terminal caps should also be specified here if capped peptides are used. Users may also specify a specific amino acid at a constrained site at the beginning of the run.  Examples: <br>`1 3 16`<br> `= 4:LEU`<br> `1:ACE 16:NME`<br> |
| Constraints | Grouped-site positional constraint amino acid pool | Amino acids that are allowed to be placed on the sites of grouped-site positional constraint. |
| Constraints | Grouped-site positional constraint | Sequence sites that are collectively assigned amino acids from the grouped-site amino acid pool without repetition. |
| Constraints | Constraints on the number of certain amino acids | Compositional constraints on number of selected amino acids. Examples: <br>`MET=0`<br> `ALA>=2`<br> `ASN<=5`<br> `PHE=1`<br> |

### Execute
Once the input files are prepared, users can submit the following example script to run PepAD:

>
	#!/bin/bash
    
    /path_to_PepAD/src/PepAD

### Analyze output files
**Table 2:** PepAD output files 
|File extension| Description|
|----------------|------------|
|`Energy profile`| Records Steps, Sequences, Γ<sub>score</sub>, ΔG<sub>binding</sub>, $$\lambda \times P_{aggregation}$$, and RMSD |
|`Energy details`| Records Steps, Trials, Sequences, Γ<sub>score</sub>, ΔG<sub>binding</sub>, ΔE<sub>binding</sub>, ΔTS<sub>conf</sub>, $$\lambda \times P_{aggregation}$$, MC moves, and Trial results|
|`Minimal energy`| Records Steps, Trials, Γ<sub>score</sub>, and Sequences with minimal Γ<sub>score</sub> along the search|
|`PDB files`     | PDB files for peptides with minimal scores during the evolution|
#### terms
- Γ<sub>score</sub>: PepAD score function, unit: kcal/mol

- ΔG<sub>binding</sub>: Binding free energy, unit: kcal/mol

- ΔE<sub>binding</sub>: Binding energy, unit: kcal/mol

- ΔTS<sub>conf</sub>: Configuration entropy, unit: kcal/mol

- P<sub>aggregation</sub>: Aggregation propensity, unit: kcal/mol

- $$\lambda$$: Aggregation propensity weighting factor, unit: 1

- RMSD: Total Root Mean Square Displacement (RMSD) of a β-sheet, unit: Angstrom.

### PepAD Analysis Code
The analysis code ranks the peptides from the lowest score, removes duplicates, and plots score versus steps and sheet RMSD versus step. This code is not required for peptide design, but it offers a convenient approach to identifying promising peptides.
#### Requirements
- Python 3
- [NumPy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [Matplotlib](https://matplotlib.org/)

A `YAML` file is included for creating a `conda` environment with the required packages.

#### Arguments
- `--top [n]` : Report *n* unique peptides with the lowest scores ranked in order.  
- `--plot [score|rmsd|both|none]` : Plot score vs step, RMSD vs step, both, or disable plots.  
- `--score_rolling [i]` : Smooth the score plot using a rolling average over *i* steps.  
- `--rmsd_rolling [i]` : Smooth the RMSD plot using a rolling average over *i* steps.  

#### Example
Once PepAD outputs are generated, user can submit the following example script to run PepAD analysis:
>
	conda activate pepad_analysis
	python pepad_analyzer.py --top 10 --plot both --score_rolling 100 --rmsd_rolling 100

output:
>
```text
--- 10 unique peptides with best score (energy profile) ---
 step  Sequence   Score   G_bind   Pagg   rmsd  Counts
  823  IFKVMGD   -19.09  -19.04   0.05   0.10       1
  802  VWKVMGE   -18.87  -18.97  -0.10   0.10       5
 3630  VIKIVGD   -18.79  -18.97  -0.18   0.17      11
 3631  VIKIMGD   -18.60  -18.86  -0.25   0.17       6
  707  KIVIVGD   -18.50  -18.79  -0.29   0.10       5
  829  VFKVEGM   -18.36  -18.53  -0.17   0.10       2
 2058  KIMIVGD   -18.36  -18.72  -0.37   0.22      12
  828  IFKVEGM   -18.30  -18.45  -0.16   0.10       3
 2083  VIKIMGE   -18.28  -18.85  -0.57   0.22      13
  711  KVIIVGD   -18.18  -18.47  -0.29   0.10       3
```

<p align="center"> 
  <img width="800" alt="step_evolution" src="https://github.com/user-attachments/assets/15786422-7758-4cfb-a3e5-2ff3c10b5856" />
  </p>
<p align="center"><b>Fig. 5. The Γ<sub>score</sub> and RMSD evolution in a 10,000-step PepAD run of designing 7-mer antiparallel peptides</b>.</figcaption></p>

### Initial Structure Builder
We provide a supplemental **initial structure builder** to build two β-sheet fibril backbones suitable for PepAD. The builder can create the β-cross spine in one of the eight steric zipper classes [4].  
#### Requirements
- Python 3
- "RotamerLibrary" in the code working directory
- [PeptideBuilder](https://github.com/mtien/PeptideBuilder) — for constructing peptides in parallel or antiparallel β-strand configuration [12]  
- [Hydride](https://github.com/biotite-dev/hydride) — for adding hydrogens at chemically reasonable positions [13]

A `YAML` file is included for creating a `conda` environment with the required packages.

#### Arguments
- `-seq`: Peptide sequence in upper case one-letter code. (required)
- `-c`: Fibril class, an integer from 1 to 8, corresponding to the 8 classes of steric zipper. (required)
- `-sh`: One β-sheet shifts along the x-axis by a space of -1, 0, or 1 residue.
- `-n`: Number of strands per sheet.
- `-p`: Terminal patch (caps). "0" = no patches, "1" =  patch with ACE and NME, "2" = patch with ACE and NHE.
- `-r`: Residue packing mode. "e" = even-numbered residues inside, "o" = odd-numbered residues inside".
- `-f`: Output PDB file format. "0" = PepAD format and "1" = AMBER format.
- `-d1`: Strand-strand distance (Å).
- `-d2`: Sheet-sheet distance (Å).
- `-d3`: One β-sheet shifts along the y-axis (fibril axis) by a distance (Å).
- `-o`: Output file name.
#### Example
User can submit the following example script to run build an initial structure:
>
	conda activate buildpep
	python Initial_structure_builder.py -seq "GNNQQNY" \
    -c "1" -sh "-1" -n "4" -p "2" -r "e" -f "1" \
    -d1 "4.8" -d2 "10" -d3 "2.4" -o "pep2"

output:
<p align="center"> 
  <img width="800" alt="step_evolution" src="https://github.com/user-attachments/assets/2c485625-af0e-4c3f-8a94-31fa2de56725" />
  </p>
<p align="center"><b>Fig. 6. A 7-mer parallel peptide backbone generated by the initial structure builder.</b>  .</figcaption></p>




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

[12]	M.Z. Tien, D.K. Sydykova, A.G. Meyer, C.O. Wilke, PeptideBuilder: A simple Python library to generate model peptides, PeerJ 1 (2013) e80. https://doi.org/10.7717/peerj.80.

[13]	P. Kunzmann, J.M. Anter, K. Hamacher, Adding hydrogen atoms to molecular models via fragment superimposition, Algorithms Mol. Biol. 17 (2022) 7. https://doi.org/10.1186/s13015-022-00215-x.


