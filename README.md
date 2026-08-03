# PepAD

Peptide Assembly Design (PepAD) is a Monte Carlo (MC)–based algorithm for discovering amyloid-forming peptides that self-assemble into user-defined reference structures [1–3]. PepAD designs peptide by optimizing the binding affinity within a predefined fibril backbone, which is composed of two stacked β-sheets.

PepAD allows users to design peptides with specified sequence length, amino acid composition, and positional or compositional constraints. During the design process, residues on the peptide backbone are modified through MC moves. The new sequence is evaluated using a physics-based score function that quantifies the peptide self-assembly propensity.

To run PepAD, users can use the precompiled container or compile the code in an
HPC environment with an Intel Fortran compiler. Users then provide an initial
fibril structure in PDB format and a text-based input file that defines parameters
and constraints. PepAD outputs peptide sequences with their corresponding scores,
as well as structural files for the top-scoring designs, which can be used in
downstream simulations or experiments.

**Current version: v1.42.** This version uses named `PARAMETER = value` entries instead of the fixed, position-dependent input format. Inputs written for v1.37 are not compatible.

## What PepAD does

PepAD starts from an initial peptide fibril structure, typically consisting of two stacked β-sheets corresponding to a specific class of cross-β spine [4,5]. A sequence is placed on the fixed peptide backbone and optimized through three Monte Carlo moves:

1. **Residue mutation** replaces a residue with another one
2. **Residue exchange** swaps two residues in the sequence.
3. **Sheet position perturbation** translates one β-sheet while the other
   β-sheet remains fixed.

<p align="center">
  <img width="450" alt="Simplified PepAD algorithm flow chart" src="https://github.com/user-attachments/assets/b0a786f8-a65d-48dd-a610-cdfaaca9e7fb" />
</p>
<p align="center"><b>Fig. 1.</b> Simplified flow chart describing the PepAD algorithm.</p>

After each move, PepAD evaluates

$$ Γ_{score} = ΔG_{binding} - \lambda \times P_{aggregation} $$

More negative scores indicate stronger predicted binding at the evaluated configuration. The binding term is calculated with an MM/GBSA-based model [6–8], and the aggregation term follows the Zyggregator approach [9–11].

## Repository layout

| Path | Contents |
| --- | --- |
| [`src/main_v1.42-2.f90`](src/main_v1.42-2.f90) | PepAD source code |
| [`src/input.example.txt`](src/input.example.txt) | Annotated input template |
| `src/lib/` | Force-field and rotamer data |
| [`examples/`](examples/) | Examples |
| [`Initial_structures/`](Initial_structures/) | Prepared fibril backbones and supporting files |
| [`Initial_structure_builder/`](Initial_structure_builder/) | Initial Structure Builder |
| [`PepAD_analyzer/`](PepAD_analyzer/) | PepAD output-analysis tools |
| [`archive/`](archive/) | Old source code, documentation, analyzer, and examples |

## Running PepAD

PepAD can be run with a precompiled container or compiled in the local
environment.

## Method 1: Use Docker image with Apptainer

The container includes a precompiled PepAD executable. Users do not need an Intel Fortran compiler for this method.

### 1. Download the SIF file

Load Apptainer if it is provided as a module on the local system. Go to the directory where the container will be stored and pull the Docker image from the GitHub Container Registry:

```bash
module load apptainer  # if required on the local system
mkdir -p /path/to/PepAD_container
cd /path/to/PepAD_container
apptainer pull PepAD_package.sif docker://ghcr.io/carolhall-ncsu-cbe/pepad-user-friendly-package:1.42-2
```

This creates `PepAD_package.sif` in `/path/to/PepAD_container`.

### 2. Prepare a run directory

Each run needs:

- an initial structure PDB file
- `input.txt`.

Then create a run directory and change to the run directory. Copy the input template into the run directory as `input.txt`, copy the initial
structure, and edit `input.txt`:

Example: if a run directory is created in the `/PepAD-User-Friendly-Package`
```bash
mkdir -p run1
cp src/input.example.txt run1/input.txt
cp Initial_structures/comp1/comp1.pdb run1/
cd run1
```

### 3. Run PepAD

Bind the run directory to `/work` inside the container. PepAD reads `input.txt` and writes its output files in the run directory:

```bash

module load apptainer  # if required on the local system
PEPAD_SIF=/path/to/PepAD_container/PepAD_package.sif  # Path to PepAD_package.sif

RUN_DIR="$(pwd -P)"    # The real path of the working directory
apptainer run --bind "${RUN_DIR}:/work" --pwd /work "$PEPAD_SIF" # execute PepAD

```
- `--bind "${RUN_DIR}:/work"` maps the current run directory on the host to `/work` inside the container.
- `--pwd /work` sets `/work` as the working directory inside the container.

### 4. Check the manual

```bash
module load apptainer  # if required on the local system
PEPAD_SIF=/path/to/PepAD_container/PepAD_package.sif

apptainer run "$PEPAD_SIF" --help
```
or
```bash
apptainer exec "$PEPAD_SIF" man PepAD
```

When using a batch script, load the Apptainer module in the script if it is required on the local system.

## Method 2: Download and compile PepAD

<p align="center">
  <img width="1200" alt="PepAD peptide-design workflow" src="https://github.com/user-attachments/assets/982d5573-7dac-4381-a9ce-651992b16658" />
</p>
<p align="center"><b>Fig. 2.</b> A workflow chart describing how to design peptides.</p>

### 1. Download and compile PepAD

Download or clone this repository. PepAD is written in Fortran 90 and requires the Intel Fortran compiler `ifx`. The provided script `compile_code.sh`, which locates the `src/`, compiles the PepAD code: 

```bash
git clone https://github.com/CarolHall-NCSU-CBE/PepAD-User-Friendly-Package.git
cd PepAD-User-Friendly-Package
bash src/compile_code.sh
```

By default, the script creates this installation:

```text
src/PepAD/
|-- PepAD
|-- lib/
`-- man/
    `-- man1/
        `-- PepAD.1
```

The PepAD executable is therefore `src/PepAD/PepAD`, not `src/PepAD`. The runtime `lib/` directory should be at same level with the PepAD executable.

To install PepAD somewhere else, provide the installation directory as the first argument:

```bash
bash src/compile_code.sh /path/to/PepAD_install
```

### 2. Add PepAD to `PATH`

For the default installation, run these commands from `PepAD-User-Friendly-Package/`:

```bash
PEPAD_INSTALL="$(pwd -P)/src/PepAD"
export PATH="$PEPAD_INSTALL:$PATH"
export MANPATH="$PEPAD_INSTALL/man:${MANPATH:-}"
```

These commands affect the current shell. To make them permanent, add the two `export` commands, using the absolute installation path, to the shell startup file.

Test whether the executable is working by referring to manual:

```bash
module load PrgEnv-intel # Load Intel Compilers if necessary
PepAD --help
man PepAD
```

### 3. Prepare a run directory

Each run needs:

- an initial structure PDB file
- `input.txt`.

Then create a run directory and change to the run directory. Copy the input template into the run directory as `input.txt`, copy the initial
structure, and edit `input.txt`:

Example: if a run directory is created in the `PepAD-User-Friendly-Package/`
```bash
mkdir -p run1
cp src/input.example.txt run1/input.txt
cp Initial_structures/comp1/comp1.pdb run1/
cd run1
```

### 4. Run PepAD

In `run1`, execute PepAD. It reads `input.txt` and writes the results in the run directory:

```bash
PepAD
```

If PepAD has not been added to `PATH`, the default installation can instead be run from `run1` using:

```bash
../src/PepAD/PepAD
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
- Comments must begin with `#` or `!`, including after a value.
- Each parameter may appear at most once.
- Unknown parameters and invalid or duplicate values stop the run with an error.
- List values may be separated by spaces or commas.
- The only required parameters are `PDBFILE` and `N_STEPS`.

A minimal input is:

```text
PDBFILE = comp1.pdb
N_STEPS = 10000
```

PepAD automatically determines the number of chains, number of amino acids per chain, energy groups, β-sheets, and composition from the PDB file. Other parameters have default values.

### Initial structure and restart settings

| Parameter | Meaning | Default |
| --- | --- | --- |
| `PDBFILE` | Initial structure filename | **Required** |
| `N_AA` | Amino acids per chain, excluding caps | Auto-recognized |
| `N_CHAINS` | Number of peptide chains | Auto-recognized |
| `GROUP_1`, `GROUP_2` | Chain IDs in the two binding-energy groups | Auto-recognized |
| `SHEET_1`, `SHEET_2`, ... | Chain IDs in each β-sheet | Auto-recognized |
| `RESTART` | `0` for a fresh run; `1` to restart | `0` |

- User-defined values for `N_AA`, `N_CHAINS`, `GROUP_1`, `GROUP_2`, and `SHEET_n` override the values automatically detected from the input PDB file.
- `N_AA` and `N_CHAINS` must either both be provided or both be omitted.
- If custom energy groups are specified, both `GROUP_1` and `GROUP_2` are required, and every PDB chain must appear exactly once across the two groups.
- Custom `SHEET_n` parameters must be consecutively numbered beginning with `SHEET_1`, and every chain must appear exactly once across all defined sheets.

Example:

```text
PDBFILE   = comp1.pdb # REQUIRED
N_AA       = 7        # Can be omitted
N_CHAINS   = 8        # Can be omitted
GROUP_1    = 1 3 5 7  # Can be omitted
GROUP_2    = 2 4 6 8  # Can be omitted
SHEET_1    = 1 2 3 4  # Can be omitted
SHEET_2    = 5 6 7 8  # Can be omitted
RESTART    = 0        # Can be omitted
```

### Monte Carlo and energy settings

| Parameter | Meaning | Default |
| --- | --- | ---: |
| `RANSEED` | Integer random seed | System clock |
| `N_STEPS` | Number of Monte Carlo steps | **Required** |
| `N_TOP` | Number of lowest-scoring unique peptide recorded | 10 |
| `KBT_SEQ` | Sequence-move kBT when annealing is off | `1.0` |
| `KBT_SHEETMOVE` | Sheet-move kBT when annealing is off | `0.6` |
| `KBT_SEQ_HIGH` | Initial sequence kBT when annealing is on | `2.0` |
| `KBT_SEQ_LOW` | Final sequence kBT when annealing is on | `0.5` |
| `KBT_SHEETMOVE_HIGH` | Initial sheet-move kBT when annealing is on | `1.2` |
| `KBT_SHEETMOVE_LOW` | Final sheet-move kBT when annealing is on | `0.3` |
| `ANNEAL_STAGES` | `0` disables annealing; a positive integer sets its stage count | `0` |
| `RMSD_MAX` | Maximum sheet-move RMSD in angstroms | `3.0` |
| `PAGG_WEIGHT` | Aggregation-propensity weight, lambda | `2.0` |
| `SHEETMOVE` | `0` disables or `1` enables sheet perturbations | `1` |
| `STEPS_BEFORE_SHEETMOVE` | Sequence steps before the first sheet move | `500` |
| `STEPS_BETWEEN_SHEETMOVE` | Sequence steps between sheet moves | `200` |
| `ESURF_MODE` | Nonpolar solvation mode: `0` off, `1` full ARVO calculation, `2` cached local ARVO calculation using a neighbor list | `0` |

- When `ANNEAL_STAGES = 0`, PepAD uses `KBT_SEQ` and `KBT_SHEETMOVE`.
- When annealing is enabled, it instead uses the corresponding `HIGH` and `LOW` values.
- All kBT values must be positive, and each high value must be greater than or equal to its low value.
- ARVO is the geometric algorithm used to calculate solvent-accessible surface area. `ESURF_MODE = 1` performs a full calculation after each MC move. `ESURF_MODE = 2` reuses cached values and recalculates only the local region affected by a mutation, making the non-polar solvation energy calculation slightly faster. In current PepAD version, both `ESURF_MODE = 1` and `ESURF_MODE = 2` are both extremely computationally expensive.

### Amino acid compositions

The input uses minimum and maximum counts to define ranges for four amino acid types:

| Parameter | Meaning | Default |
| --- | --- | --- |
| `N_HYDROPHOBIC_MIN`, `N_HYDROPHOBIC_MAX` | Minimum and maximum number of hydrophobic amino acids: `GLY ALA VAL LEU ILE MET PHE TYR TRP` | Composition of the input structure |
| `N_POLAR_MIN`, `N_POLAR_MAX` | Minimum and maximum number of polar amino acids: `SER ASN GLN THR HIE` | Composition of the input structure |
| `N_CHARGED_MIN`, `N_CHARGED_MAX` | Minimum and maximum number of charged amino acids: `ARG LYS GLU ASP` | Composition of the input structure |
| `N_OTHER_MIN`, `N_OTHER_MAX` | Minimum and maximum number of other amino acids: `CYS PRO` | Composition of the input structure |

- If all eight parameters are omitted, PepAD uses the same amino-acid composition as the initial structure.
- When designing compositions different from the initial structure, each category must specify a minimum, a maximum, or both.
- An omitted minimum is set to `0`.
- An omitted maximum is set to `N_AA`.
- Equal minimum and maximum values impose an exact count.
- Number of amino acids in each types include residues in single-site positional constraints and grouped-site positional constraints, and exclude `ACE`, `NME`, and `NHE` caps.

Example 1: Designing 7-mer peptides with dynamic amino acid compositions

```text
N_HYDROPHOBIC_MIN = 0
N_HYDROPHOBIC_MAX = 7
N_POLAR_MIN       = 0
N_POLAR_MAX       = 7
N_CHARGED_MIN     = 0
N_CHARGED_MAX     = 2
N_OTHER_MIN       = 0
N_OTHER_MAX       = 0
```

Example 2: Designing 7-mer peptides with fixed compositions:

```text
N_HYDROPHOBIC_MIN = 5
N_HYDROPHOBIC_MAX = 5
N_POLAR_MIN       = 0
N_POLAR_MAX       = 0
N_CHARGED_MIN     = 2
N_CHARGED_MAX     = 2
N_OTHER_MIN       = 0
N_OTHER_MAX       = 0
```

### Compositional constraints

Any amino acid can also have its own count range:

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

- Enter only a minimum, only a maximum, or both.
- Equal MIN and MAX values impose an exact count.
- Amino acids without these parameters remain unconstrained, subject to their category limits.
- \<AA\> can be replaced with `GLY`, `ALA`, `VAL`, `LEU`, `ILE`, `MET`, `PHE`, `TYR`, `TRP`, `SER`, `ASN`, `GLN`, `THR`, `HIE`, `ARG`, `LYS`, `GLU`, `ASP`, `CYS`, and `PRO`.
<br>

### Positional constraints

- Cap positions (`ACE`, `NME`, and `NHE`) are recognized automatically and should not be listed as constraints. 
- Site numbers here refer only to amino-acid positions.

#### Single-site positional constraint

`SINGLE_SITE_CONSTRAINTS` accepts following values separated by spaces:

- `site_num` keeps the amino acid already present at that site.
- `site_num:AA` first changes the site to the named amino acid, then keeps it fixed.

**Do not place spaces around the colon:**

Example:
```text
SINGLE_SITE_CONSTRAINTS = 5:ASN 6
```

Use `NONE` or omit the parameter when no single-site positional constraint is needed.

#### Grouped-site positional constraint

Grouped-site constraints assign amino acids to selected positions from a user-defined pool. Example:

```text
GROUPED_SITE_AA_POOL     = ASN SER THR ALA VAL
GROUPED_SITE_CONSTRAINTS = 1 2 3
```
Here, sites 1, 2, and 3 are assigned amino acids from the `GROUPED_SITE_AA_POOL`, and each amino acid can be selected up to once.

PepAD allows duplicate amino acids in the pool. Example: 
```text
GROUPED_SITE_AA_POOL     = ASN SER ALA ALA VAL
GROUPED_SITE_CONSTRAINTS = 2 3 5
```
Here, `ALA` may appear on sites 2, 3, and 5 at most twice.

- The number of amino acids in `GROUPED_SITE_AA_POOL` must not be less than the number of sites in `GROUPED_SITE_CONSTRAINTS`.
- To allow an amino acid to appear one additional time, enter it one additional time in the pool.
- `GROUPED_SITE_CONSTRAINTS` cannot overlap with sites listed in `SINGLE_SITE_CONSTRAINTS`.
- PepAD checks whether the grouped-site pool is compatible with the overall composition limits and compositional constraints.


### Annotated input.txt template

The template is [`src/input.example.txt`](src/input.example.txt). Use it as the starting point for new designs.

## Initial structure requirements

Users must provide an initial structure PDB file for PepAD. The structure usually consists of individual peptides in a configuration of cross-β spine. Prepared structures are available in [`Initial_structures/`](Initial_structures/). Two recommended approaches for preparing the initial structure are as follows.

First, users can search for existing amyloid fibril structures in the [RCSB Protein Data Bank](https://www.rcsb.org/). For example, the Sup35 prion segment GNNQQNY structure (PDB ID: 2OMM) can be replicated to generate a larger Class-1 cross-β spine. The resulting structure should be relaxed using explicit-solvent atomistic molecular dynamics simulation to remove atomic overlaps and then converted to PepAD's required PDB format.

<p align="center">
  <img width="650" alt="Preparing a PepAD structure from a Protein Data Bank fibril" src="https://github.com/user-attachments/assets/c476b3f8-90dd-411e-837f-1fc97d476770" />
</p>
<p align="center"><b>Fig. 3.</b> Prepare initial structure using existing fibril structure from Protein Data Bank.</p>

Second, users can build artificial amyloid backbones using molecular modeling tools such as UCSF Chimera and Packmol. For example, Aβ(16–22) KLVFFAE was built with UCSF Chimera, packed into a two-layer Class-8 antiparallel β-sheet fibril using Packmol, and relaxed by explicit-solvent molecular dynamics simulation. To reduce PepAD computational cost, the middle eight peptides were extracted from the relaxed fibril and saved as `comp1.pdb`. The amyloid backbone can also be generated using the Initial Structure Builder provided in this package.

<p align="center">
  <img width="800" alt="Preparing an artificial amyloid backbone with molecular modeling tools" src="https://github.com/user-attachments/assets/3b1b6026-a947-4a42-98e6-0e8600c831d6" />
</p>
<p align="center"><b>Fig. 4.</b> Prepare initial structure using molecular modeling tool.</p>

## Output files

PepAD writes output into the run directory.

| Output | Contents |
| --- | --- |
| `energyprofile.txt` | Steps, Sequences, Γ<sub>score</sub>, ΔG<sub>binding</sub>, $$\lambda \times P_{aggregation}$$, and RMSD |
| `energydetails.txt` | Steps, Trials, Sequences, Γ<sub>score</sub>, ΔE<sub>VDW</sub>, ΔE<sub>ELE</sub>,ΔE<sub>SGB</sub>,ΔE<sub>SURF</sub>, ΔTS<sub>conf</sub>, $$\lambda P_{aggregation}$$, ΔΔE<sub>VDW</sub>, ΔΔE<sub>ELE</sub>, ΔΔE<sub>SGB</sub>, ΔΔE<sub>SURF</sub>, ΔΔTS<sub>conf</sub>, Δ $$\lambda P_{aggregation}$$, MC moves, and Trial results|
| `minimum_energy.txt` | Lowest-score sequences encountered during the search |
| `backup4backbone.txt` | Restart information |
| `pdbfiles/` | Structures saved for minimum-score designs |

Run settings and progress are printed to standard output; an HPC submission script can redirect that stream to a job report. Energy quantities are reported in kcal/mol. RMSD is reported in angstroms.

## PepAD tools
We provide a PepAD-related Python package that contains two modules: `builder` and `analyzer`. These tools are helpful but are not required to run PepAD.
- `builder`: create the β-cross spine in one of the ten steric zipper classes [4,13].
- `analyzer`: analyze PepAD results (energy profile and energy details) and derive best-scoring peptides, plotting score change over steps, as well as how individual energy terms contribute to the MC acceptance or rejection.

### Install PepAD tools
The installation only needs to be done once.
Clone the PepAD repository if it has not already been downloaded:
 ```bash
git clone https://github.com/CarolHall-NCSU-CBE/PepAD-User-Friendly-Package.git
cd PepAD-User-Friendly-Package
 ```

Create and activate a Conda environment, then install PepAD tools via `pip`:
```bash
conda create -n pepad-tools python=3.12 pip
conda activate pepad-tools
python -m pip install ./PepAD_tools/
```
On HPC systems, users may create the environment to a specified path
```bash
conda create -p /path/to/env_pepad_tools python=3.12 pip
conda activate /path/to/env_pepad_tools
python -m pip install ./PepAD_tools/
```

Both builder and the analyzer can be called in the "pepad-tools" environment
```bash
builder --help
analyzer --help
```
### Load PepAD Tools Later
After opening a new terminal session, users need to activate the environment again:
```bash
conda activate pepad-tools
```
On HPC systems
```bash
conda activate /path/to/env_pepad_tools
```

### Use Initial Structure Builder (builder)
#### Arguments
- `-h`, `--help`: Show help message.
- `-s SEQUENCE`, `--seq SEQUENCE`, `--seqA SEQUENCE`: Input peptide sequence. (required)
- `-c {1,2,3,4,5,6,7,8,9,10}`, `--class {1,2,3,4,5,6,7,8,9,10}`: Integer from 1 to 10, corresponding to the 10 classes of steric zipper.
- `-dx DX`, `--dx DX`: Strand-strand distance (Å) along x. (default: 4.8 Å)
- `-dz DZ`, `--dz DZ`: Sheet-sheet distance (Å) along z. (default: 11.5 Å)
- `-x X`, `--x X`: Sheet-2 shifts in the x direction in units of `0.5 * dx`. (default: 0.0)
- `-y Y`, `--y Y`: Sheet-2 shifts in the y direction in units of residue spacing (3.465 Å). (default: 0.0)
- `-n N_CHAINS`, `--chains N_CHAINS`: Integer. Strands per sheet. (default: 8)
- `-r {e,o}`, `--core {e,o}`: Packing mode. `e` = even-numbered residues packed inside; `o` = odd-numbered residues packed inside. (default: e)
- `-g {minus,plus}`, `--registry {minus,plus}`: Relative direction of registry shift inside antiparallel beta-sheets. (default: minus)
- `-f {0,1}`, `--format {0,1}`: PDB format. `0` = PepAD format; `1` = AMBER format. (default: 0)
- `-C {0,1,2}`, `--cap {0,1,2}`: Terminal caps. `0` = uncapped; `1` = ACE+NME; `2` = ACE+NHE. (default: 0)
- `-o OUTPUT`, `--output OUTPUT`: Output PDB filename. (default: peptide_sheets)

#### Example
User can submit the following example script to run build an initial structure:
```bash
builder --seq GNNQQNY --class 1 --dx 4.8 --dz 11.0 --x -1 --y 1 --chains 4 --cap 2 --core e --format 1 --output pep2.pdb
```
output:
<p align="center">
  <img width="800" alt="Parallel 7-mer peptide backbone generated by the Initial Structure Builder" src="https://github.com/user-attachments/assets/2c485625-af0e-4c3f-8a94-31fa2de56725" />
</p>
<p align="center"><b>Fig. 5.</b> A 7-mer parallel peptide backbone generated by the Initial Structure Builder.</p>

### Use PepAD Analyzer (analyzer)
#### Arguments
- `-h`, `--help`: Show help message.
- `--directory DIRECTORY`: PepAD run directory. (default: current directory)
- `--profiles`: Analyze `energyprofile.txt` and write `PepAD_report.txt`.
- `--top TOP`: Integer. Number of top-scoring unique peptides written to `PepAD report.txt`. (default: 10)
- `--details`: Analyze `energydetails.txt` and write `Detail_report.txt`.
- `--plot {score,rmsd,gbind,pagg} [...]`: Generate one or more optional step-wise plots. Requires `--profiles`.
- `--violin {energy2, contrib2, energy5, contrib5} [...]`: Generate violin plots showing energy-term contributions to the negative part of the score change (Δscore) in accepted trials or to the positive part of the (Δscore) in rejected trials. Requires `--details`.
  - `energy2`: Plot the distributions of ΔΔG_binding and Δ(-P_agg). (default)
  - `contrib2`: Plot their percentage contributions to the negative part of Δscore in accepted trials or to the positive part of Δscore in rejected trials.
  - `energy5`: Plot the distributions of ΔΔE_VDW, ΔΔ(E_ELE + G_GB), ΔΔG_SURF, Δ(-TS), and Δ(-P_agg).
  - `contrib5`: Plot the corresponding percentage contributions of these five terms.
- `--rolling ROLLING`: Integer. Number of rolling-average windows for plots. (default: 0)

At least one of `--profiles` or `--details` is needed. Both can be used in the same command.

#### Example
Once PepAD outputs are generated, user can submit the following example script to run PepAD analysis:
>
	analyzer --top 10 --plot score rmsd --rolling 100

output:
PepAD_report.txt
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
Score and RMSD vs steps plot
<p align="center">
  <img width="800" alt="PepAD score and RMSD evolution" src="https://github.com/user-attachments/assets/15786422-7758-4cfb-a3e5-2ff3c10b5856" />
</p>
<p align="center"><b>Fig. 6.</b> The Γ<sub>score</sub> and RMSD evolution in a 10,000-step PepAD run of designing 7-mer antiparallel peptides.</p>

## References

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
13. Stroud, J. C. The Zipper Groups of the Amyloid State of Proteins. Acta Cryst D 2013, 69 (4), 540–545. https://doi.org/10.1107/S0907444912050548.

    
## License

This project is distributed under the terms in [`LICENSE`](LICENSE).

## Acknowledgment

This project is supported by the National Science Foundation under Award No.
1931430, “Element: Computational Toolkit to Discover Peptides that Self-assemble
into User-selected Structures,” and Award No. 2347712, “CDS&E: Computational
Design of Peptide-Based Biorecognition Elements.” Any opinions, findings, and
conclusions or recommendations expressed in this project are those of the authors
and do not necessarily reflect the views of the National Science Foundation.
