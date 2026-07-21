# PepAD

Peptide Assembly Design (PepAD) is a Monte Carlo sequence-design program for
discovering peptides predicted to self-assemble on a user-supplied amyloid fibril
backbone. PepAD explores peptide sequences and sheet arrangements, evaluates each
trial with a physics-based score, and reports low-scoring designs for downstream
simulation or experimental study.

**Current version: v1.42.** The input format changed in v1.42 from a fixed,
position-dependent form to named `FLAG = value` entries. Old v1.37 inputs are not
compatible with v1.42. The v1.37 source, documentation, analyzer, and examples are
preserved in [`archive/v1.37`](archive/v1.37).

## What PepAD does

PepAD begins with a fibril model, usually two stacked beta sheets in one of the
cross-beta steric-zipper classes. A sequence is placed on the fixed peptide
backbone and optimized through three Monte Carlo moves:

1. **Residue mutation** replaces a residue with another residue from the same
   hydration category.
2. **Residue exchange** swaps two positions in the sequence.
3. **Sheet position perturbation** translates one beta sheet while the other
   sheet remains fixed.

After each move, PepAD evaluates

$$
\Gamma_{\mathrm{score}} = \Delta G_{\mathrm{binding}}
- \lambda P_{\mathrm{aggregation}}.
$$

More-negative scores indicate stronger predicted binding at the evaluated
configuration. The binding term is calculated with an MM/GBSA-based model, and
the aggregation term follows the Zyggregator approach.

## Repository layout

| Path | Contents |
| --- | --- |
| [`src/main.f90`](src/main.f90) | Current PepAD v1.42 source |
| [`src/input.example.txt`](src/input.example.txt) | Annotated v1.42 input template |
| `src/lib/` | Runtime force-field and rotamer data |
| [`examples/`](examples/) | Location for v1.42 examples |
| [`Initial structures/`](Initial%20structures/) | Prepared fibril backbones and supporting files |
| [`Initial structure builder/`](Initial%20structure%20builder/) | Optional Python tool for building fibril backbones |
| [`archive/v1.37/`](archive/v1.37/) | Complete legacy v1.37 material |

The repository root tracks the current release. Version-specific historical
files belong in `archive/`; do not place old source files beside `src/main.f90`.

## Quick start

### 1. Compile PepAD

PepAD is written in Fortran 90 and is intended for an HPC environment with an
Intel Fortran compiler. From `src/`, compile the program beside the `lib/`
directory:

```bash
module load PrgEnv-intel
cd src
ifx -O2 -o PepAD main.f90
```

Older Intel environments may use `ifort` instead of `ifx`. The executable finds
its parameter library relative to its own location, so keep `PepAD` in `src/`
beside `lib/` unless both are deployed together.

To list all supported v1.42 input flags and defaults:

```bash
./PepAD --help
```

### 2. Prepare a run directory

Each run needs:

- an initial fibril structure in PDB format; and
- a text file named exactly `input.txt`.

Copy [`src/input.example.txt`](src/input.example.txt) into the run directory,
rename it to `input.txt`, and edit it for the selected PDB structure:

```bash
cd /path/to/PepAD-User-Friendly-Package
mkdir my_run
cp src/input.example.txt my_run/input.txt
cp "Initial structures/comp1/comp1.pdb" my_run/
cd my_run
```

### 3. Run PepAD

Launch the executable from the run directory so that PepAD reads that directory's
`input.txt` and writes its results there:

```bash
../src/PepAD
```

Adjust the executable path when the run directory is elsewhere.

## v1.42 input format

### Syntax rules

Every nonblank input line has the form:

```text
FLAG = value
```

- Flag order does not matter.
- Flag names are case-insensitive; uppercase is recommended.
- Blank lines are allowed.
- Comments may begin with `#` or `!`, including after a value.
- Each flag may appear at most once.
- Unknown flags and invalid or duplicate values stop the run with an error.
- List values may be separated by spaces; chain lists also accept commas.
- The only required flags are `PDBFILE` and `N_STEPS`.

A minimal input is:

```text
PDBFILE = comp1.pdb
N_STEPS = 10000
```

With this minimal form, PepAD determines chain length, chain count, energy
groups, beta sheets, and composition from the PDB file. Other calculations use
the defaults listed below.

### Initial structure and restart settings

| Flag | Meaning | Default or rule |
| --- | --- | --- |
| `PDBFILE` | Initial structure filename | **Required**; normally relative to the run directory |
| `N_RESIDUES` | Residues per chain, including caps | Auto-detected when omitted |
| `N_CHAINS` | Total peptide chains | Auto-detected when omitted |
| `GROUP_1`, `GROUP_2` | Chain IDs in the two binding-energy groups | DSSP-based assignment when both are omitted |
| `SHEET_1`, `SHEET_2`, ... | Chain IDs in each beta sheet | DSSP-based assignment when all are omitted |
| `RESTART` | `0` for a fresh run; `1` to restart | `0` |

`N_RESIDUES` and `N_CHAINS` must either both be present or both be absent. If
custom energy groups are provided, both `GROUP_1` and `GROUP_2` are required and
every PDB chain must occur exactly once across the pair. Custom `SHEET_n` flags
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

| Flag | Meaning | Default |
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

### Amino-acid composition ranges

The v1.42 input uses minimum and maximum counts instead of the four exact counts
used by v1.37:

| Category | Amino acids | Flags |
| --- | --- | --- |
| Hydrophobic | `GLY ALA VAL LEU ILE MET PHE TYR TRP` | `N_HYDROPHOBIC_MIN`, `N_HYDROPHOBIC_MAX` |
| Polar | `SER ASN GLN THR HIE` | `N_POLAR_MIN`, `N_POLAR_MAX` |
| Charged | `ARG LYS GLU ASP` | `N_CHARGED_MIN`, `N_CHARGED_MAX` |
| Other | `CYS PRO` | `N_OTHER_MIN`, `N_OTHER_MAX` |

If all eight category flags are omitted, PepAD uses the starting PDB composition.
If custom composition mode is used, every category must have a minimum, a
maximum, or both. An omitted minimum becomes `0`; an omitted maximum becomes
`N_RESIDUES`. Equal minimum and maximum values impose an exact count.

Composition totals exclude `ACE`, `NME`, and `NHE` caps, but include residues at
single-site and grouped-site constrained positions.

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

### Specific amino-acid count constraints

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
count. Amino acids without these flags remain unconstrained, subject to their
category limits. Supported names are `GLY`, `ALA`, `VAL`, `LEU`, `ILE`, `MET`,
`PHE`, `TYR`, `TRP`, `SER`, `ASN`, `GLN`, `THR`, `HIE`, `ARG`, `LYS`, `GLU`,
`ASP`, `CYS`, and `PRO`.

### Positional constraints

Cap positions (`ACE`, `NME`, and `NHE`) are recognized automatically and should
not be listed as constraints. Site numbers refer only to amino-acid positions;
caps are skipped when numbering. Thus, amino acid 3 is chain site 4 when chain
site 1 is an ACE cap.

#### Single-site constraints

`SINGLE_SITE_CONSTRAINTS` accepts space-separated entries:

- `site` keeps the amino acid already present at that site.
- `site:AA` first changes the site to the named amino acid, then keeps it fixed.

Do not place spaces around the colon:

```text
SINGLE_SITE_CONSTRAINTS = 5:ASN 6
```

Use `NONE` or omit the flag when no single-site constraint is needed.

#### Grouped-site constraints

Grouped sites draw residues from a user-defined pool:

```text
GROUPED_SITE_AA_POOL     = ASN SER THR ALA VAL
GROUPED_SITE_CONSTRAINTS = 1 2 3
```

The pool must contain at least as many entries as there are grouped sites. A
residue may be repeated in the pool to make more copies available; each remaining
copy has equal selection probability. Grouped sites may not overlap caps or
single-site constraints. PepAD validates the pool against the category and
specific-amino-acid count limits.

### Complete annotated template

The maintained template is [`src/input.example.txt`](src/input.example.txt). Use
it as the starting point for new calculations. Do not adapt a v1.37 input by
changing only its labels: the old positional format and exact composition fields
have different semantics.

## Initial structure requirements

The PDB file should contain peptide chains arranged as a cross-beta fibril and use
atom and residue naming compatible with the parameter files in `src/lib/`.
Prepared structures are available in [`Initial structures/`](Initial%20structures/).
An initial structure may be obtained by:

1. adapting an experimental amyloid structure from the
   [RCSB Protein Data Bank](https://www.rcsb.org/) and relaxing it to remove
   clashes; or
2. building a steric-zipper model with the supplied
   [`Initial structure builder`](Initial%20structure%20builder/), followed by an
   appropriate relaxation procedure.

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

## Legacy v1.37 calculations

For exact reproduction of an older run, use the source, input template, analyzer,
and examples together from [`archive/v1.37`](archive/v1.37). Do not combine a
v1.37 input with the v1.42 executable. For long-term release management, create a
Git tag or GitHub Release for each stable version; the default branch should
continue to contain only the current version.

## Citation and references

If PepAD contributes to published work, cite the relevant PepAD design papers:

1. Sarma, S. et al. Design of parallel beta-sheet nanofibrils using Monte Carlo
   search, coarse-grained simulations, and experimental testing. *Protein
   Science* **33** (2024), e5102. https://doi.org/10.1002/pro.5102
2. Xiao, X. et al. Sequence patterns and signatures: Computational and
   experimental discovery of amyloid-forming peptides. *PNAS Nexus* **1**
   (2022), pgac263. https://doi.org/10.1093/pnasnexus/pgac263
3. Xiao, X. et al. De novo design of peptides that coassemble into beta
   sheet-based nanofibrils. *Science Advances* **7** (2021), eabf7668.
   https://doi.org/10.1126/sciadv.abf7668

## License

This project is distributed under the terms in [`LICENSE`](LICENSE).
