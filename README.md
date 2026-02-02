# PepAD-User-Friendly-Package
## Introduction  
 
Peptide Assembly Design (PepAD) is a Monte Carlo (MC)–based algorithm for discovering amyloid-forming peptides that self-assemble into user-defined reference structures [1–3]. PepAD designs peptide by optimizing the binding affinity within a predefined fibril backbone, which is composed of two stacked β-sheets.

PepAD allows users to design peptides with specified sequence length, amino acid composition, and positional or compositional constraints. During the design process, residues on the peptide backbone are modified through MC moves. The new sequence is evaluated using a physics-based score function that quantifies the peptide self-assembly propensity.

To run PepAD, users first need to compile the code in HPC with Intel@ Fortran compiler. Then need to provide an initial fibril structure (PDB format) and a text-based input file that defines parameters and constraints. PepAD outputs peptide sequences with their corresponding scores, as well as structural files for the top-scoring designs, which can be further used in downstream simulations or experiments.


## PepAD algorithm
  PepAD starts from an initial peptide fibril structure, typically consisting of two stacked β-sheets corresponding to a specific class of cross-β spine[4-5]. At beign, a random sequence is draped upon the backbone (–NH–Cα–CO–), and its binding affinity is evaluated using a score function (Γ<sub>score</sub>). PepAD then explores new sequences using three MC moves—residue mutation, residue exchange, and sheet position perturbation. After each move, the Γ<sub>score</sub> is evaluated for the new configuration, and the new design is accepted or rejected based on Metropolis criterion. This process is iterated over multiple steps, and peptides with the lowest Γ<sub>score</sub> are identified as the best designs.

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
<p align="center">
  <img width="1778" height="491" alt="image" src="https://github.com/user-attachments/assets/9fe5fec4-7abc-4b22-81e6-e7443452672b" />
  </p>
<p align="center"><b>Fig. 3.</b> A workflow chart describing how to design peptides.</figcaption></p>


### Compile code
The coding language is Fortran90. Users need to load Intel Fortran Compiler (ifx) for compilation. In the command line interface, navigate to the `/src/`, and use the following command:
 
>

	# Example compilation using the Intel Fortran compiler on a Linux-based system
	
	# Users should modify these commands according to their local environment. 

	module load PrgEnv-intel         # load Intel compiler

    ifx -o PepAD main.f90            #  compile PepAD code


> **Note:** Code must be compiled in the same directory as `/lib/`, since parameter files are required at runtime. After compilation, the PepAD can be run from any directory.

### Set up input files
PepAD requires an `input.txt` file and an `initial_structure.pdb` to start searching. 

#### initial structure
Users must provide an initial structure PDB file for PepAD. The PDB file format must match that in `/example/`. The structure usually consists of inidividual peptides in a configuration of cross-β spine. Two recommended approaches for preparing the initial structure are as follows: 

(a) Search for existing amyloid fibril structures in the Protein Data Bank ([RCSB PDB: Homepage](https://www.rcsb.org/)). Users need to modify the format of the PDB file to match PepAD's required input format (see Fig. 4, top).

(b) Build artificial amyloid backbones using the peptide-building tool provided in this package, or prepare with other molecular modelling tools (UCSF Chimera or Packmol) (see Fig. 4, bottom).

<p align="center"> 
  <img width="450" alt="image" src="https://github.com/user-attachments/assets/e6c6a3dd-9f5e-4f17-b85f-c81f36759425" />
  <img width="450" alt="image" src="https://github.com/user-attachments/assets/1311ca29-03ca-4b8c-8014-9b195bbd81c7" />

  </p>
<p align="center"><b>Fig. 4. (a) comp-1.pdb: A 7-mer fibril structure built using crystal structure of peptide GNNQQNY (PDB ID: 2omm) and molecular dynamics simulation . (b) comp-2.pdb: A 14-mer fibril strucutre built using peptide-building tool in PepAD package.</b></figcaption></p>

#### input.txt
The 'input.txt' contains the necessary parameters for designing a fibril structure.
**Table 1:** Input paremeters for PepAD
| Item | Description |
|------|-------------|
| Number of residues per chain (including caps) | An integer number. Total number of residues, including terminal residues. |
| Number of chains in total | An integer number. Total number of peptide chains. |
| Input PDB file name | A string. Name of the initial structure PDB file. |
| Recalculation | Switch for restarting PepAD runs:<br>`0` # fresh run<br>`1` # restart from last completed step |
| Group 1 | Indices of peptide chains in group 1.|
| Group 2 | Indices of peptide chains in group 2.|
| Number of MC moves | An integer number. Total steps for performing MC moves. |
| Ekt(sequence) | A real number. The temperature factor kBT is used in the Metropolis criterion for residue mutation/exchange. |
| Ekt(sheet move) | A real number. The temperature factor kBT is used in the Metropolis criterion for the sheet position perturbation move. |
| Allowed maximum RMSD (ang) | A real number. The maximum displacement of one β-sheet can shift relative to another β-sheet. |
| Aggregation propensity weight | A real number. The weighting factor λ adjusts the contribution of aggregation propensity in the score function. |
| Hydrophobic | Number of hydrophobic residues (G, A, L, V, I, M, F, Y, W). |
| Polar | Number of polar residues (N, Q, S, T). |
| Charged | Number of charged residues (R, K, H, E, D). |
| Other | Number of residues that are not classified as hydrophobic, polar, or charged, (C, P). |
| Single-site positional constraint | The sites on which amino acids do not change. |
| Grouped-site positional constraint amino acids pool | A set of amino acids that could only be placed on the grouped positional constraint sites. |
| Grouped-site positional constraint | The sites on which amino acids can only be selected from the set of amino acids defined previously, without repetition. |
| Constraints on the number of certain amino acids | Limit the count of specific amino acids. For examples:<br>`MET=0` # No Met<br>`Ala>=2` # The peptide must have at least 2 Ala<br>`Asn<=5` # The peptide must have no more than 5 Asn |

### Submit a job
Once the input files are prepared, users can submit the following example script to run PepAD:

>
	#!/bin/bash

    mkdir -p pdbfiles
    
    /path_to_PepAD/src/PepAD

### Output files
**Table 2:** PepAD output files 
|File extension| Description|
|----------------|------------|
|`Energy profile`| Records Steps, Sequences, Γ<sub>score</sub>, ΔG<sub>bind</sub>, $$\lambda \times P_{agg}$$, and RMSD |
|`Energy details`| Records Steps, Trials, Sequences, Γ<sub>score</sub>, ΔG<sub>bind</sub>, ΔE<sub>bind</sub>, ΔTS<sub>conf</sub>, $$\lambda \times P_{agg}$$, MC moves, and Trial results|
|`Minimal energy`| Records Steps, Trials, Γ<sub>score</sub>, and Sequences with minimal Γ<sub>score</sub> along the search|
|`PDB files`     | PDB files for peptides with minimal scores during the evolution|
#### terms
- Γ<sub>score</sub>: PepAD score function

- ΔG<sub>binding</sub>: Binding free energy

- ΔE<sub>binding</sub>: Binding energy

- ΔTS<sub>conf</sub>: Configuration entropy

- P<sub>aggregation: Aggregation propensity

- $$\lambda$$: Aggregation propensity weighting factor

- RMSD: Total Root Mean Square Displacement (RMSD) of a β-sheet.

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
	---10 unique peptides with best score (energy profile)---
	 step Sequence  Score  E_bind  S_bind  I_hydrophobic  I_propensity  G_bind  Pagg  rmsd_x  rmsd_y  rmsd_z  rmsd  Counts
	 3917  VWKVVGD -20.17  -22.94   -3.06          -1.29          1.44  -19.88  0.15    0.20    0.00    0.74  0.76      22
	 3911  VWKVVGE -20.03  -23.51   -3.45          -1.11          1.10  -20.06 -0.01    0.20    0.00    0.74  0.76      30
	 6816  VWKVMGD -19.97  -23.96   -4.21          -1.31          1.42  -19.75  0.11    0.37    0.09    0.83  0.92      43
	 4006  KWVVVGD -19.73  -23.41   -3.86          -1.29          1.38  -19.55  0.09    0.20    0.00    0.60  0.63      17
	 6836  VWKVMGE -19.71  -24.34   -4.53          -1.13          1.08  -19.81 -0.05    0.37    0.09    0.83  0.92      23
	 1887  KWMVVGD -19.38  -23.28   -4.00          -1.31          1.36  -19.27  0.05    0.01    0.00    0.80  0.80      22
	 5064  KWMVVGE -19.33  -24.09   -4.55          -1.13          1.02  -19.54 -0.11    0.33    0.06    0.68  0.76       8
	 1895  KWVVMGD -19.10  -23.82   -4.82          -1.31          1.36  -19.00  0.05    0.01    0.00    0.80  0.80      30
	 3918  KWVVVGE -19.03  -23.74   -4.57          -1.11          1.04  -19.17 -0.07    0.20    0.00    0.74  0.76      15
	 1131  KWIIDGG -18.86  -22.63   -3.65          -1.23          1.17  -18.98 -0.06    0.02    0.00    0.64  0.64      24


<p align="center"> 
  <img width="800" alt="step_evolution" src="https://github.com/user-attachments/assets/9f3abf04-6600-43bc-95ac-8cab8fa2f0a0" />
  </p>
<p align="center"><b>Fig. 5. The Γ<sub>score</sub> and RMSD evolution in a 10,000-step PepAD run of designing 7-mer antiparallel peptides.</b>  .</figcaption></p>

### Initial Structure Builder
We provide a supplemental **initial structure builder** to build two β-sheet fibril backbones suitable for PepAD. The builder can create the β-cross spine in one of the eight steric zipper classes defined by Sawaya and Eisenberg [4].  
#### Requirements
- Python 3 
- [Biopython](https://biopython.org/) — for storing peptide structures [12]  
- [PeptideBuilder](https://github.com/mtien/PeptideBuilder) — for constructing peptides in parallel or antiparallel β-strand configuration [13]  
- [Hydride](https://github.com/biotite-dev/hydride) — for adding hydrogens at chemically reasonable positions [14]

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

[12]	P.J.A. Cock, T. Antao, J.T. Chang, B.A. Chapman, C.J. Cox, A. Dalke, I. Friedberg, T. Hamelryck, F. Kauff, B. Wilczynski, M.J.L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics 25 (2009) 1422–1423. https://doi.org/10.1093/bioinformatics/btp163.

[13]	M.Z. Tien, D.K. Sydykova, A.G. Meyer, C.O. Wilke, PeptideBuilder: A simple Python library to generate model peptides, PeerJ 1 (2013) e80. https://doi.org/10.7717/peerj.80.

[14]	P. Kunzmann, J.M. Anter, K. Hamacher, Adding hydrogen atoms to molecular models via fragment superimposition, Algorithms Mol. Biol. 17 (2022) 7. https://doi.org/10.1186/s13015-022-00215-x.


