# PepAD Tools

PepAD Tools provides two terminal programs:

- `builder` creates initial beta-sheet peptide structures.
- `analyzer` reads PepAD simulation output and generates reports and plots.

The builder supports uniform steric-zipper classes 1-10, hybrid types 1-6,
terminal caps, and one- or two-sequence sheets.

## Install With Conda

Clone the PepAD repository if it has not already been downloaded:

```bash
git clone https://github.com/CarolHall-NCSU-CBE/PepAD-User-Friendly-Package.git
cd PepAD-User-Friendly-Package
```

Create an environment and install PepAD Tools:

```bash
conda create -n pepad-tools python=3.12 pip
conda activate pepad-tools
python -m pip install ./PepAD_tools/
```

The installation only needs to be performed once. Activate the environment
again whenever a new terminal or batch job is started:

```bash
conda activate pepad-tools
```

## Install Without Conda

Linux or macOS:

```bash
cd PepAD-User-Friendly-Package
python3 -m venv ~/env_pepad
source ~/env_pepad/bin/activate
python -m pip install ./PepAD_tools/
```

Windows:

```powershell
cd PepAD-User-Friendly-Package
py -3.12 -m venv env_pepad
env_pepad\Scripts\activate
python -m pip install ./PepAD_tools/
```

NumPy, pandas, Matplotlib, and the rotamer templates are installed with the
package.

## Builder

Show all terminal inputs:

```bash
builder --help
```

Build a Class 1 sheet:

```bash
builder --seq GNNQQNY --class 1 --chains 4 --output peptide.pdb
```

### Builder Arguments

- `-s`, `--seq`: Peptide sequence in one-letter amino-acid code. This is
  sequence A when `--seqB` is supplied. Required.
- `-b`, `--seqB`: Optional second peptide sequence. It must have the same
  length as sequence A.
- `-p1`, `--pattern1`: Sheet-1 A/B pattern, such as `ABAB`. Required with
  `--seqB`.
- `-p2`, `--pattern2`: Sheet-2 A/B pattern. Required with `--seqB`.
- `-c`, `--class`: Uniform steric-zipper class from 1 to 10.
- `-H`, `--hybrid`: Parallel-antiparallel hybrid type from 1 to 6. Supply
  either `--class` or `--hybrid`, but not both.
- `-dx`, `--dx`: Strand-strand distance along x in angstroms. Default: `4.8`.
- `-dz`, `--dz`: Sheet-sheet distance along z in angstroms. Default: `11.5`.
- `-x`, `--x`: Sheet-2 shift along x in half-strand-spacing units. For
  example, `--x 1` shifts sheet 2 by `dx/2`.
- `-y`, `--y`: Sheet-2 shift along y, the fibril axis, in residue-spacing
  units. One unit is 3.465 angstroms.
- `-n`, `--chains`: Number of strands per sheet for a one-sequence build.
  Default: `8`. Pattern lengths determine strand counts for two sequences.
- `-r`, `--core`: Residue packing mode. `e` packs even-numbered residues
  inside; `o` packs odd-numbered residues inside.
- `-g`, `--registry`: Registry direction, `minus` or `plus`, for
  antiparallel sheets.
- `-f`, `--format`: PDB format. `0` is PepAD format; `1` is AMBER format.
- `-C`, `--cap`: Terminal caps. `0` is uncapped, `1` is ACE+NME, and `2` is
  ACE+NHE.
- `-o`, `--output`: Output PDB filename. The `.pdb` suffix is added when
  absent.
- `-t`, `--tilt`: Maximum angle between the middle-residue CA-to-CB vector
  and the z-axis for parallel-angle strands. Default: `10` degrees.
- `--no-tilt`: Disable the middle-residue CA-to-CB angle limit.

### Builder Example

Build a four-strand-per-sheet Class 1 structure with ACE and NHE caps:

```bash
conda activate pepad-tools
builder --seq GNNQQNY --class 1 --x -1 --chains 4 \
    --cap 2 --core e --format 1 --dx 4.8 --dz 10 --output pep2
```

The command can be run from any working directory. Generated PDB files are
written to that directory.

The same interface is also available as:

```bash
python -m pepad_tools --help
```

## Analyzer

The analyzer reads parameters from `input.txt` and initial structure PDB
file. It analyzes `energyprofile.txt`, `energydetails.txt`, or both. Run it
from the PepAD run directory, or select that directory with
`--directory`.

Show all analyzer inputs:

```bash
analyzer --help
```

Write the best 100 unique sequences from `energyprofile.txt`:

```bash
analyzer --profiles --top 100
```

Write `Detail_report.txt` without figures:

```bash
analyzer --details
```

Write both reports and generate every optional plot:

```bash
analyzer --profiles --top 100 \
    --plot score rmsd gbind pagg \
    --details --violin energy contribution
```

## Test

From the repository folder after installation:

```bash
python -m unittest discover -s tests -v
```
