# PepAD

Peptide Assembly Design (PepAD) is a Monte Carlo sequence-design program for
discovering peptides predicted to self-assemble on a user-supplied amyloid fibril
backbone. PepAD explores peptide sequences and sheet arrangements, evaluates each
trial with a physics-based score, and reports low-scoring designs for downstream
simulation or experimental study.

**Current version: v1.42.** This version uses named `PARAMETER = value` entries
instead of the fixed, position-dependent input format. Inputs written for v1.37
are not compatible.

## What PepAD does

PepAD starts from an initial peptide fibril structure, typically consisting of
two stacked β-sheets corresponding to a specific class of cross-β spine [4,5]. A
sequence is placed on the fixed peptide backbone and optimized through three
Monte Carlo moves:

1. **Residue mutation** replaces a residue with another residue from the same
   hydration category.
2. **Residue exchange** swaps two positions in the sequence.
3. **Sheet position perturbation** translates one beta sheet while the other
   sheet remains fixed.

<p align="center">
  <img width="450" alt="Simplified PepAD algorithm flow chart" src="https://github.com/user-attachments/assets/b0a786f8-a65d-48dd-a610-cdfaaca9e7fb" />
</p>
<p align="center"><b>Fig. 1.</b> Simplified flow chart describing the PepAD algorithm.</p>

After each move, PepAD evaluates

$$ Γ_{score} = ΔG_{binding} - \lambda \times P_{aggregation} $$

More-negative scores indicate stronger predicted binding at the evaluated
configuration. The binding term is calculated with an MM/GBSA-based model
[6–8], and the aggregation term follows the Zyggregator approach [9–11].

## Repository layout

| Path | Contents |
| --- | --- |
| [`src/main_v1.42.f90`](src/main_v1.42.f90) | PepAD source code |
| [`src/input.example.txt`](src/input.example.txt) | Annotated input template |
| `src/lib/` | Runtime force-field and rotamer data |
| [`examples/`](examples/) | Examples |
| [`Initial structures/`](Initial%20structures/) | Prepared fibril backbones and supporting files |
| [`Initial structure builder/`](Initial%20structure%20builder/) | Initial Structure Builder |
| [`archive/v1.37/`](archive/v1.37/) | Old source code, documentation, analyzer, and examples |

The old files are in [`archive/v1.37/`](archive/v1.37/).

## Quick start

<p align="center">
  <img width="1200" alt="PepAD peptide-design workflow" src="https://github.com/user-attachments/assets/982d5573-7dac-4381-a9ce-651992b16658" />
</p>
<p align="center"><b>Fig. 2.</b> A workflow chart describing how to design peptides.</p>

### 1. Compile PepAD

PepAD is written in Fortran 90 and is intended for an HPC environment with an
Intel Fortran compiler. From `src/`, compile the program beside the `lib/`
directory:

```bash
module load PrgEnv-intel
cd /path/to/PepAD-User-Friendly-Package/src
ifx -O2 -o PepAD main_v1.42.f90
```

Older Intel environments may use `ifort` instead of `ifx`. The executable finds
its parameter library relative to its own location, so keep `PepAD` in `src/`
beside `lib/` unless both are deployed together.

To list all supported input parameters and defaults:

```bash
/path/to/PepAD-User-Friendly-Package/src/PepAD --help
```

### 2. Prepare a run directory

Each run needs:

- an initial fibril structure in PDB format; and
- a text file named exactly `input.txt`.

Create a working directory outside the PepAD package. Copy
[`src/input.example.txt`](src/input.example.txt) into the working directory as
`input.txt`, copy the initial PDB structure, and edit `input.txt`:

```bash
mkdir -p /path/to/PepAD_runs/my_run
cp /path/to/PepAD-User-Friendly-Package/src/input.example.txt /path/to/PepAD_runs/my_run/input.txt
cp "/path/to/PepAD-User-Friendly-Package/Initial structures/comp1/comp1.pdb" /path/to/PepAD_runs/my_run/
cd /path/to/PepAD_runs/my_run
```

### 3. Run PepAD

Remember the path to the compiled PepAD executable. From the working directory,
run PepAD using its full path. PepAD reads `input.txt` and writes the results in
the working directory:

```bash
/path/to/PepAD-User-Friendly-Package/src/PepAD
```

## Input format

### Syntax rules

Every nonblank input line has the form:

```text
PARAMETER = value
```

- Parameter order does not matter.
- Parameter names are case-insensitive; uppercase is recommended.
- Blank lines are allowed.
- Comments may begin with `#` or `!`, including after a value.
- Each parameter may appear at most once.
- Unknown parameters and invalid or duplicate values stop the run with an error.
- List values may be separated by spaces; chain lists also accept commas.
- The only required parameters are `PDBFILE` and `N_STEPS`.

A minimal input is:

```text
PDBFILE = comp1.pdb
N_STEPS = 10000
```

With this minimal form, PepAD determines chain length, chain count, energy
groups, beta sheets, and composition from the PDB file. Other calculations use
the defaults listed below.

### Initial structure and restart settings

| Parameter | Meaning | Default or rule |
| --- | --- | --- |
| `PDBFILE` | Initial structure filename | **Required**; normally relative to the run directory |
| `N_RESIDUES` | Residues per chain, including caps | Auto-detected when omitted |
| `N_CHAINS` | Total peptide chains | Auto-detected when omitted |
| `GROUP_1`, `GROUP_2` | Chain IDs in the two binding-energy groups | DSSP-based assignment when both are omitted |
| `SHEET_1`, `SHEET_2`, ... | Chain IDs in each beta sheet | DSSP-based assignment when all are omitted |
| `RESTART` | `0` for a fresh run; `1` to restart | `0` |

`N_RESIDUES` and `N_CHAINS` must either both be present or both be absent. If
custom energy groups are provided, both `GROUP_1` and `GROUP_2` are required and
every PDB chain must occur exactly once across the pair. Custom `SHEET_n` parameters
must be consecutively numbered from `SHEET_1`, and every chain must occur exactly
once across all sheets.

Example:

```text
PDBFILE   = comp1.pdb
N_RESIDUES = 7
N_CHAINS   = 8
GROUP_1    = 1 3 5 7
GROUP_2    = 2 4 6 8
SHEET_1    = 1 2 3 4
SHEET_2    = 5 6 7 8
RESTART    = 0
```

### Monte Carlo and energy settings

| Parameter | Meaning | Default |
| --- | --- | ---: |
| `RANSEED` | Integer random seed | System clock |
| `N_STEPS` | Number of Monte Carlo steps | **Required** |
| `KBT_SEQ` | Sequence-move kBT when annealing is off | `1.0` |
| `KBT_SHEETMOVE` | Sheet-move kBT when annealing is off | `0.6` |
| `KBT_SEQ_HIGH`, `KBT_SEQ_LOW` | Initial and final sequence kBT during annealing | `2.0`, `0.5` |
| `KBT_SHEETMOVE_HIGH`, `KBT_SHEETMOVE_LOW` | Initial and final sheet-move kBT during annealing | `1.2`, `0.3` |
| `ANNEAL_STAGES` | `0` disables annealing; a positive integer sets its stage count | `0` |
| `RMSD_MAX` | Maximum sheet-move RMSD in angstroms | `3.0` |
| `PAGG_WEIGHT` | Aggregation-propensity weight, lambda | `2.0` |
| `SHEETMOVE` | `0` disables or `1` enables sheet perturbations | `1` |
| `STEPS_BEFORE_SHEETMOVE` | Sequence steps before the first sheet move | `500` |
| `STEPS_BETWEEN_SHEETMOVE` | Sequence steps between sheet moves | `200` |
| `ESURF_MODE` | Nonpolar solvation mode: `0` off, `1` full ARVO, `2` cached/incremental ARVO | `0` |

When `ANNEAL_STAGES = 0`, PepAD uses `KBT_SEQ` and `KBT_SHEETMOVE`. When
annealing is enabled, it instead uses the corresponding `HIGH` and `LOW` values.
All kBT values must be positive, and each high value must be greater than or
equal to its low value.

### Amino acid compositions

The input uses minimum and maximum counts instead of the four exact counts
used by v1.37:

| Category | Amino acids | Parameters |
| --- | --- | --- |
| Hydrophobic | `GLY ALA VAL LEU ILE MET PHE TYR TRP` | `N_HYDROPHOBIC_MIN`, `N_HYDROPHOBIC_MAX` |
| Polar | `SER ASN GLN THR HIE` | `N_POLAR_MIN`, `N_POLAR_MAX` |
| Charged | `ARG LYS GLU ASP` | `N_CHARGED_MIN`, `N_CHARGED_MAX` |
| Other | `CYS PRO` | `N_OTHER_MIN`, `N_OTHER_MAX` |

If all eight category parameters are omitted, PepAD uses the starting PDB composition.
If custom composition mode is used, every category must have a minimum, a
maximum, or both. An omitted minimum becomes `0`; an omitted maximum becomes
`N_RESIDUES`. Equal minimum and maximum values impose an exact count.

Composition totals exclude `ACE`, `NME`, and `NHE` caps, but include residues at
sites defined by single-site positional constraints and grouped-site positional
constraints.

Example allowing a range of polar and charged content while fixing the
hydrophobic count:

```text
N_HYDROPHOBIC_MIN = 4
N_HYDROPHOBIC_MAX = 4
N_POLAR_MIN       = 0
N_POLAR_MAX       = 2
N_CHARGED_MIN     = 0
N_CHARGED_MAX     = 2
N_OTHER_MIN       = 0
N_OTHER_MAX       = 0
```

### Compositional constraints

Any standard PepAD amino acid can also have its own count range:

```text
N_<AA>_MIN = integer
N_<AA>_MAX = integer
```

For example:

```text
N_ASN_MIN = 0
N_ASN_MAX = 1
N_SER_MIN = 1
N_SER_MAX = 3
```

Supply only a minimum, only a maximum, or both. Equal values impose an exact
count. Amino acids without these parameters remain unconstrained, subject to their
category limits. Supported names are `GLY`, `ALA`, `VAL`, `LEU`, `ILE`, `MET`,
`PHE`, `TYR`, `TRP`, `SER`, `ASN`, `GLN`, `THR`, `HIE`, `ARG`, `LYS`, `GLU`,
`ASP`, `CYS`, and `PRO`.

### Positional constraints

Cap positions (`ACE`, `NME`, and `NHE`) are recognized automatically and should
not be listed as constraints. Site numbers refer only to amino-acid positions;
caps are skipped when numbering. Thus, amino acid 3 is chain site 4 when chain
site 1 is an ACE cap.

#### Single-site positional constraint

`SINGLE_SITE_CONSTRAINTS` accepts space-separated entries:

- `site` keeps the amino acid already present at that site.
- `site:AA` first changes the site to the named amino acid, then keeps it fixed.

Do not place spaces around the colon:

```text
SINGLE_SITE_CONSTRAINTS = 5:ASN 6
```

Use `NONE` or omit the parameter when no single-site positional constraint is needed.

#### Grouped-site positional constraint

Grouped sites draw residues from a user-defined pool:

```text
GROUPED_SITE_AA_POOL     = ASN SER THR ALA VAL
GROUPED_SITE_CONSTRAINTS = 1 2 3
```

The pool must contain at least as many entries as there are grouped sites. A
residue may be repeated in the pool to make more copies available; each remaining
copy has equal selection probability. Grouped sites may not overlap caps or
single-site positional constraints. PepAD validates the pool against the amino
acid compositions and compositional constraints.

### Complete annotated template

The maintained template is [`src/input.example.txt`](src/input.example.txt). Use
it as the starting point for new calculations. Do not adapt a v1.37 input by
changing only its labels: the old positional format and exact composition fields
have different semantics.

## Initial structure requirements

Users must provide an initial structure PDB file for PepAD. The structure usually
consists of individual peptides in a configuration of cross-β spine. Prepared
structures are available in [`Initial structures/`](Initial%20structures/). Two
recommended approaches for preparing the initial structure are as follows.

First, users can search for existing amyloid fibril structures in the
[RCSB Protein Data Bank](https://www.rcsb.org/). For example, the Sup35 prion
segment GNNQQNY structure (PDB ID: 2OMM) can be replicated to generate a larger
Class-1 cross-β spine. The resulting structure should be relaxed using
explicit-solvent atomistic molecular dynamics simulation to remove atomic
overlaps and then converted to PepAD's required PDB format.

<p align="center">
  <img width="650" alt="Preparing a PepAD structure from a Protein Data Bank fibril" src="https://github.com/user-attachments/assets/c476b3f8-90dd-411e-837f-1fc97d476770" />
</p>
<p align="center"><b>Fig. 3.</b> Prepare initial structure using existing fibril structure from Protein Data Bank.</p>

Second, users can build artificial amyloid backbones using molecular modeling
tools such as UCSF Chimera and Packmol. For example, Aβ(16–22) KLVFFAE was built
with UCSF Chimera, packed into a two-layer Class-8 antiparallel β-sheet fibril
using Packmol, and relaxed by explicit-solvent molecular dynamics simulation. To
reduce PepAD computational cost, the middle eight peptides were extracted from
the relaxed fibril and saved as `comp1.pdb`. The amyloid backbone can also be
generated using the Initial Structure Builder provided in this package.

<p align="center">
  <img width="800" alt="Preparing an artificial amyloid backbone with molecular modeling tools" src="https://github.com/user-attachments/assets/3b1b6026-a947-4a42-98e6-0e8600c831d6" />
</p>
<p align="center"><b>Fig. 4.</b> Prepare initial structure using molecular modeling tool.</p>

## Output files

PepAD writes output into the run directory.

| Output | Contents |
| --- | --- |
| `energyprofile.txt` | Step, sequence, score, binding free energy, weighted aggregation propensity, and sheet RMSD |
| `energydetails.txt` | Detailed trial energies, Monte Carlo move, and acceptance result |
| `minimum_energy.txt` | Lowest-score sequences encountered during the search |
| `backup4backbone.txt` | Restart information |
| `pdbfiles/` | Structures saved for minimum-score designs |

Run settings and progress are printed to standard output; an HPC submission
script can redirect that stream to a job report. Energy quantities are reported
in kcal/mol. RMSD is reported in angstroms.

<p align="center">
  <img width="800" alt="PepAD score and RMSD evolution" src="https://github.com/user-attachments/assets/15786422-7758-4cfb-a3e5-2ff3c10b5856" />
</p>
<p align="center"><b>Fig. 5.</b> The Γ<sub>score</sub> and RMSD evolution in a 10,000-step PepAD run of designing 7-mer antiparallel peptides.</p>

### Initial Structure Builder

We provide a supplemental **Initial Structure Builder** to build two β-sheet
fibril backbones suitable for PepAD. The Initial Structure Builder can create the
β-cross spine in one of the eight steric zipper classes [4]. It uses
PeptideBuilder to construct peptides and Hydride to add hydrogens [12,13].

<p align="center">
  <img width="800" alt="Parallel 7-mer peptide backbone generated by the Initial Structure Builder" src="https://github.com/user-attachments/assets/2c485625-af0e-4c3f-8a94-31fa2de56725" />
</p>
<p align="center"><b>Fig. 6.</b> A 7-mer parallel peptide backbone generated by the Initial Structure Builder.</p>

## Citation and references

If PepAD contributes to published work, please cite the relevant PepAD design
papers [1–3].

1. S. Sarma, T.R. Sudarshan, V. Nguyen, A.S. Robang, X. Xiao, J.V. Le,
   M.E. Helmicki, A.K. Paravastu, and C.K. Hall, “Design of parallel β-sheet
   nanofibrils using Monte Carlo search, coarse-grained simulations, and
   experimental testing,” *Protein Science* **33** (2024), e5102.
   https://doi.org/10.1002/pro.5102
2. X. Xiao, A.S. Robang, S. Sarma, J.V. Le, M.E. Helmicki, M.J. Lambert,
   R. Guerrero-Ferreira, J. Arboleda-Echavarria, A.K. Paravastu, and C.K. Hall,
   “Sequence patterns and signatures: Computational and experimental discovery
   of amyloid-forming peptides,” *PNAS Nexus* **1** (2022), pgac263.
   https://doi.org/10.1093/pnasnexus/pgac263
3. X. Xiao, Y. Wang, D.T. Seroski, K.M. Wong, R. Liu, A.K. Paravastu,
   G.A. Hudalla, and C.K. Hall, “De novo design of peptides that coassemble into
   β-sheet-based nanofibrils,” *Science Advances* **7** (2021), eabf7668.
   https://doi.org/10.1126/sciadv.abf7668
4. M.R. Sawaya, S. Sambashivan, R. Nelson, M.I. Ivanova, S.A. Sievers,
   M.I. Apostol, M.J. Thompson, M. Balbirnie, J.J.W. Wiltzius, H.T. McFarlane,
   A.Ø. Madsen, C. Riekel, and D. Eisenberg, “Atomic structures of amyloid
   cross-β spines reveal varied steric zippers,” *Nature* **447** (2007),
   453–457. https://doi.org/10.1038/nature05695
5. R. Nelson, M.R. Sawaya, M. Balbirnie, A.Ø. Madsen, C. Riekel, R. Grothe,
   and D. Eisenberg, “Structure of the cross-β spine of amyloid-like fibrils,”
   *Nature* **435** (2005), 773–778. https://doi.org/10.1038/nature03680
6. H. Gohlke, C. Kiel, and D.A. Case, “Insights into protein–protein binding by
   binding free energy calculation and free energy decomposition for the
   Ras–Raf and Ras–RalGDS complexes,” *Journal of Molecular Biology* **330**
   (2003), 891–913. https://doi.org/10.1016/S0022-2836(03)00610-7
7. G. Rastelli, A.D. Rio, G. Degliesposti, and M. Sgobba, “Fast and accurate
   predictions of binding free energies using MM-PBSA and MM-GBSA,” *Journal of
   Computational Chemistry* **31** (2010), 797–810.
   https://doi.org/10.1002/jcc.21372
8. X. Xiao, B. Zhao, P.F. Agris, and C.K. Hall, “Simulation study of the ability
   of a computationally designed peptide to recognize target tRNA Lys3 and other
   decoy tRNAs,” *Protein Science* **25** (2016), 2243–2255.
   https://doi.org/10.1002/pro.3056
9. G.G. Tartaglia, A.P. Pawar, S. Campioni, C.M. Dobson, F. Chiti, and
   M. Vendruscolo, “Prediction of aggregation-prone regions in structured
   proteins,” *Journal of Molecular Biology* **380** (2008), 425–436.
   https://doi.org/10.1016/j.jmb.2008.05.013
10. G.G. Tartaglia and M. Vendruscolo, “The Zyggregator method for predicting
    protein aggregation propensities,” *Chemical Society Reviews* **37** (2008),
    1395–1401. https://doi.org/10.1039/B706784B
11. A.P. Pawar, K.F. DuBay, J. Zurdo, F. Chiti, M. Vendruscolo, and C.M. Dobson,
    “Prediction of aggregation-prone and aggregation-susceptible regions in
    proteins associated with neurodegenerative diseases,” *Journal of Molecular
    Biology* **350** (2005), 379–392.
    https://doi.org/10.1016/j.jmb.2005.04.016
12. M.Z. Tien, D.K. Sydykova, A.G. Meyer, and C.O. Wilke, “PeptideBuilder: A
    simple Python library to generate model peptides,” *PeerJ* **1** (2013), e80.
    https://doi.org/10.7717/peerj.80
13. P. Kunzmann, J.M. Anter, and K. Hamacher, “Adding hydrogen atoms to molecular
    models via fragment superimposition,” *Algorithms for Molecular Biology*
    **17** (2022), 7. https://doi.org/10.1186/s13015-022-00215-x

## License

This project is distributed under the terms in [`LICENSE`](LICENSE).

## Acknowledgment

This project is supported by the National Science Foundation under Award No.
1931430, “Element: Computational Toolkit to Discover Peptides that Self-assemble
into User-selected Structures,” and Award No. 2347712, “CDS&E: Computational
Design of Peptide-Based Biorecognition Elements.” Any opinions, findings, and
conclusions or recommendations expressed in this project are those of the authors
and do not necessarily reflect the views of the National Science Foundation.
