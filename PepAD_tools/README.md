# PepAD Tools

PepAD Tools provides two terminal programs:

- `builder` creates initial beta-sheet peptide structures.
- `analyzer` reads PepAD simulation output and generates reports and plots.

The builder supports uniform steric-zipper classes 1-10, hybrid types 1-6,
terminal caps, and one- or two-sequence sheets.

## Install With Conda

From the repository's `pip` folder:

```bash
cd pip
conda create -n pepad python=3.12 pip
conda activate pepad
python -m pip install .
```

## Install Without Conda

Linux or macOS:

```bash
cd pip
python3 -m venv ~/env_pepad
source ~/env_pepad/bin/activate
python -m pip install .
```

Windows:

```powershell
cd pip
py -3.12 -m venv env_pepad
env_pepad\Scripts\activate
python -m pip install .
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
