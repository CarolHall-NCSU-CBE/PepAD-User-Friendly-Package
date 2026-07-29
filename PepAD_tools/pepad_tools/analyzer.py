#!/usr/bin/env python
# coding: utf-8

# In[114]:


import os
import pandas as pd
from matplotlib import pyplot as plt
import argparse
import math
import warnings


# In[ ]:


AMINO_ACID_CODES = set("ACDEFGHIKLMNPQRSTVWY")
AMINO_ACID_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIE", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}
CAP_NAMES = {"ACE", "NME", "NHE"}


def _validate_peptide_sequence(
    sequence: str, expected_length: int | None, location: str
) -> str:
    """
    Function:
        Validate one peptide sequence and return its normalized value.

    Parameters
    ----------
    sequence : str
        Peptide sequence using one-letter amino-acid codes.

    expected_length : int or None
        Required sequence length. None disables the length check.

    location : str
        File and line description used in an error message.

    Returns
    -------
    str
        Uppercase validated peptide sequence.
    """
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError(f"{location}: peptide sequence is missing.")

    invalid_codes = sorted(set(sequence) - AMINO_ACID_CODES)
    if invalid_codes:
        invalid_text = ", ".join(invalid_codes)
        raise ValueError(
            f"{location}: invalid amino-acid code(s): {invalid_text}."
        )

    if expected_length is not None and len(sequence) != expected_length:
        raise ValueError(
            f"{location}: peptide sequence has {len(sequence)} residues; "
            f"expected {expected_length}."
        )

    return sequence


def _validate_residue_sequence(
    residues: list[str], expected_length: int | None, location: str
) -> list[str]:
    """
    Function:
        Validate one energy-detail residue sequence.

    Parameters
    ----------
    residues : list of str
        One- or three-letter residue names, including caps.

    expected_length : int or None
        Required number of residues including caps.

    location : str
        File and line description used in an error message.

    Returns
    -------
    list of str
        Uppercase validated residue names.
    """
    residues = [residue.strip().upper() for residue in residues]
    if not residues:
        raise ValueError(f"{location}: peptide sequence is missing.")

    invalid_residues = []
    for residue in residues:
        base_name = residue[1:] if len(residue) == 4 else residue
        terminal_name = len(residue) == 4 and residue[0] in {"N", "C"}
        valid_name = residue in AMINO_ACID_NAMES or residue in CAP_NAMES
        valid_code = len(residue) == 1 and residue in AMINO_ACID_CODES
        if not valid_name and not valid_code:
            if not terminal_name or base_name not in AMINO_ACID_NAMES:
                invalid_residues.append(residue)

    if invalid_residues:
        invalid_text = ", ".join(sorted(set(invalid_residues)))
        raise ValueError(f"{location}: invalid residue name(s): {invalid_text}.")

    if expected_length is not None and len(residues) != expected_length:
        raise ValueError(
            f"{location}: peptide sequence has {len(residues)} residues; "
            f"expected {expected_length}."
        )

    return residues


def _normalize_step_numbers(
    steps: list[int], source_name: str, attempts: list[int] | None = None
) -> list[int]:
    """
    Function:
        Renumber duplicate, missing, or out-of-order step blocks.

    Parameters
    ----------
    steps : list of int
        Raw step number for every profile row or energy-detail trial.

    source_name : str
        Input filename displayed in correction messages.

    attempts : list of int or None, default=None
        Trial attempts used to identify energy-detail step blocks.

    Returns
    -------
    list of int
        Sequential logical step number for every input row.
    """
    if not steps:
        return []
    if attempts is not None and len(attempts) != len(steps):
        raise ValueError("Step and attempt arrays must have equal lengths.")

    block_starts = [0]
    if attempts is None:
        block_starts.extend(range(1, len(steps)))
    else:
        for index in range(1, len(steps)):
            new_step = steps[index] != steps[index - 1]
            reset_attempt = attempts[index] <= attempts[index - 1]
            if new_step or reset_attempt:
                block_starts.append(index)

    raw_blocks = [int(steps[index]) for index in block_starts]
    corrected_blocks = [raw_blocks[0]]
    messages = []
    for index in range(1, len(raw_blocks)):
        previous_raw = raw_blocks[index - 1]
        raw_step = raw_blocks[index]
        corrected_step = corrected_blocks[-1] + 1
        corrected_blocks.append(corrected_step)

        if raw_step == previous_raw:
            messages.append(
                f"Duplicate step {raw_step} found. This occurrence and all "
                "later steps were shifted back by 1 (step numbers +1)."
            )
        elif raw_step > previous_raw + 1:
            missing_start = previous_raw + 1
            missing_end = raw_step - 1
            missing_count = missing_end - missing_start + 1
            missing_label = (
                str(missing_start)
                if missing_count == 1
                else f"{missing_start}-{missing_end}"
            )
            messages.append(
                f"Missing step(s) {missing_label} found before step {raw_step}. "
                f"Step {raw_step} and all later steps were shifted forward "
                f"by {missing_count} (step numbers -{missing_count})."
            )
        elif raw_step < previous_raw:
            messages.append(
                f"Out-of-order step {raw_step} found after step {previous_raw}. "
                f"It was renumbered as logical step {corrected_step}."
            )

    normalized = [0] * len(steps)
    for block_index, start in enumerate(block_starts):
        end = (
            block_starts[block_index + 1]
            if block_index + 1 < len(block_starts)
            else len(steps)
        )
        block_length = end - start
        normalized[start:end] = [corrected_blocks[block_index]] * block_length

    if messages:
        print(f"Step inconsistencies found in {source_name}:")
        for message in messages:
            print(f"  {message}")
        print(f"  Corrected final step: {corrected_blocks[-1]}")

    return normalized


def analyze_pdb(pdb_file: str) -> tuple[int, int, int]:
    """
    Function:
        Determine chain length and number of chains from a PepAD PDB file.

    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file.

    Returns
    -------
    tuple of int
        Number of amino acids per chain, residues per chain including caps,
        and peptide chains.
    """
    res_names = []
    previous_resid = None

    # ACE, NME, and NHE are caps, not amino acids.
    cap_names = {"ACE", "NME", "NHE"}

    with open(pdb_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(
                ("ATOM", "HETATM")
            ):  # Ignore records that do not contain atom coordinates.
                continue

            res_name = line[17:21].strip().upper()
            resid = (line[21:22], line[22:26], line[26:27])

            if resid != previous_resid:
                res_names.append(res_name)
                previous_resid = resid

    N_cap = 0
    C_cap = 0
    chain_lengths = []
    local_resid = 0
    for res_name in res_names:
        if res_name == "ACE" or (
            len(res_name) == 4 and res_name.startswith("N")
        ):  # chain starts either "ACE" or "N<AA>"
            if res_name == "ACE":
                N_cap = 1
            local_resid = 1
            continue

        local_resid += 1
        if (
            res_name == "NME"
            or res_name == "NHE"
            or (len(res_name) == 4 and res_name.startswith("C"))
        ):  # chain ends "NME", "NHE" or "C<AA>"
            if res_name == "NME" or res_name == "NHE":
                C_cap = 1
            chain_lengths.append(local_resid)
            local_resid = 0

    n_residues = chain_lengths[0]  # number of residues
    n_aa = n_residues - N_cap - C_cap  # number of amino acids
    n_chains = len(chain_lengths)

    return n_aa, n_residues, n_chains


def read_input(basepath: str = ".") -> dict[str, object] | None:
    """
    Function:
        Read PepAD input parameters and detect peptide dimensions.

    Parameters
    ----------
    basepath : str, default='.'
        Directory containing input.txt and the input PDB file.

    Returns
    -------
    dict or None
        PepAD parameters, or None when input.txt is not found.
    """
    basepath = os.path.abspath(basepath)
    input_file = os.path.join(basepath, "input.txt")
    print("PepAD Analyzer working directory:", os.getcwd())
    if not os.path.isfile(input_file):
        print("input.txt is not found:", input_file)
        return None
    # Available parameters in v1.42 input.txt
    # INPUTS      : If values are not specified, use the defaults
    # INTEGER_FLAGS : Parameters that are input as intergers
    # FLOAT_FLAGS   : Parameters that are input as reals
    # AA3           : 3-letter amino acids
    INPUTS = {
        "RESTART": 0,
        "N_AA": -1,
        "N_CHAINS": -1,
        "N_RESIDUES": -1,
        "RANSEED": None,
        "KBT_SEQ": 1.0,
        "KBT_SEQ_HIGH": 2.0,
        "KBT_SEQ_LOW": 0.5,
        "RMSD_MAX": 3.0,
        "PAGG_WEIGHT": 2.0,
        "ANNEAL_STAGES": 0,
        "CONROT_FLAG": 0,
        "CONROT_ANGLE_LIMIT": 2.0,
        "CONROT_STEP": 0.1,
        "CONROT_CLOSURE_TOL": 0.10,
        "CONROT_RMSD_LIMIT": 0.5,
    }

    INTEGER_FLAGS = {
        "N_HYDROPHOBIC_MIN",
        "N_HYDROPHOBIC_MAX",
        "N_POLAR_MIN",
        "N_POLAR_MAX",
        "N_CHARGED_MIN",
        "N_CHARGED_MAX",
        "N_OTHER_MIN",
        "N_OTHER_MAX",
        "N_AA",
        "RESTART",
        "RANSEED",
        "N_STEPS",
        "ANNEAL_STAGES",
        "CONROT_FLAG",
    }

    FLOAT_FLAGS = {
        "KBT_SEQ",
        "KBT_SEQ_HIGH",
        "KBT_SEQ_LOW",
        "RMSD_MAX",
        "PAGG_WEIGHT",
        "CONROT_ANGLE_LIMIT",
        "CONROT_STEP",
        "CONROT_CLOSURE_TOL",
        "CONROT_RMSD_LIMIT",
    }

    AA3 = {
        "GLY",
        "ALA",
        "VAL",
        "LEU",
        "ILE",
        "MET",
        "PHE",
        "TYR",
        "TRP",
        "SER",
        "ASN",
        "GLN",
        "THR",
        "HIE",
        "HIS",
        "ARG",
        "LYS",
        "GLU",
        "ASP",
        "CYS",
        "PRO",
    }

    # Read v1.42 parameters
    parameters = INPUTS.copy()
    with open(input_file, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, start=1):  # Iterate all lines in input.txt
            for marker in ("#", "!"):
                line = line.split(marker, 1)[
                    0
                ]  # Remove inline comments beginning with # or !

            line = line.strip()  # Remove spaces at beginning and end of a string

            if not line:
                continue  # Skip blank and comment-only lines

            if "=" not in line:  # Raise error if no "=" in line
                raise ValueError(f"Line {i}: expected 'PARAMETER = value'.")

            key, value = line.split("=", 1)
            key = key.strip().upper()
            value = value.strip()

            if not value:  # Raise error if no value after "="
                raise ValueError(f"Line {i}: no value provided for {key}.")

            if key in INTEGER_FLAGS:
                parameters[key] = int(value)

            elif key in FLOAT_FLAGS:
                parameters[key] = float(value)

            elif key == "PDBFILE":
                # Remove optional quotation marks
                parameters[key] = value.strip("\"'")

    if "PDBFILE" not in parameters or "N_STEPS" not in parameters:
        raise ValueError("PDBFILE and N_STEPS are required in input.txt.")

    pdb_file = parameters["PDBFILE"]
    if not os.path.isabs(pdb_file):
        pdb_file = os.path.join(basepath, pdb_file)
    parameters["PDBFILE"] = os.path.abspath(pdb_file)
    parameters["BASE_DIR"] = basepath

    pdb_n_aa, pdb_n_residues, pdb_n_chains = analyze_pdb(
        parameters["PDBFILE"]
    )
    if parameters["N_AA"] not in (-1, pdb_n_aa):
        warnings.warn(
            f"N_AA={parameters['N_AA']} in input.txt does not match "
            f"N_AA={pdb_n_aa} detected from the PDB file. "
            "The PDB-detected value will be used.",
            UserWarning,
            stacklevel=2,
        )

    parameters["N_AA"] = pdb_n_aa
    parameters["N_RESIDUES"] = pdb_n_residues
    parameters["N_CHAINS"] = pdb_n_chains

    return parameters


# In[116]:


def read_energy_profile(
    basepath: str, expected_sequence_length: int | None = None
) -> pd.DataFrame:
    """
    Function:
        Read energyprofile.txt and remove malformed records with a warning.

    Parameters
    ----------
    basepath : str
        Directory containing energyprofile.txt.

    expected_sequence_length : int or None, default=None
        Expected peptide length. The first record defines it when omitted.

    Returns
    -------
    pandas.DataFrame
        Energy profile with step, sequence, energy, and RMSD columns.
    """
    energy_file = os.path.join(basepath, "energyprofile.txt")
    headers = ["step", "Sequence", "Score", "G_bind", "Pagg", "rmsd"]
    if not os.path.isfile(energy_file):
        raise FileNotFoundError(f"energyprofile.txt is not found: {energy_file}")

    with open(energy_file, "r", encoding="utf-8", errors="replace") as f:
        data_lines = [
            (line_number, line.split())
            for line_number, line in enumerate(f, start=1)
            if line.strip()
        ]

    if not data_lines:
        raise ValueError("Energy profile is incomplete. No data were found.")

    sequence_length = expected_sequence_length
    valid_rows = []
    removed_records = []
    for line_number, values in data_lines:
        step = None
        if values:
            try:
                step = int(values[0])
            except ValueError:
                pass
        step_label = str(step) if step is not None else "unknown"

        if len(values) != len(headers):
            reason = f"expected 6 values, found {len(values)}"
            removed_records.append((line_number, step_label, reason))
            continue
        if step is None:
            removed_records.append(
                (line_number, step_label, "invalid or missing step number")
            )
            continue

        try:
            sequence = _validate_peptide_sequence(
                values[1], sequence_length, "peptide sequence"
            )
        except ValueError as error:
            removed_records.append((line_number, step_label, str(error)))
            continue

        numbers = []
        invalid_columns = []
        for column, value in zip(headers[2:], values[2:]):
            try:
                number = float(value)
            except ValueError:
                number = float("nan")
            if not math.isfinite(number):
                invalid_columns.append(column)
            numbers.append(number)

        if invalid_columns:
            invalid_text = ", ".join(invalid_columns)
            reason = f"invalid numeric value(s) in {invalid_text}"
            removed_records.append((line_number, step_label, reason))
            continue

        if sequence_length is None:
            sequence_length = len(sequence)
        valid_rows.append([step, sequence, *numbers])

    if removed_records:
        details = "\n".join(
            f"  file line {line}, step {step}: {reason}"
            for line, step, reason in removed_records
        )
        warnings.warn(
            "Incomplete or invalid energy-profile record(s) removed:\n"
            f"{details}\n"
            "Removed records are treated as missing steps.",
            UserWarning,
            stacklevel=2,
        )

    if not valid_rows:
        raise ValueError("Energy profile contains no complete records.")

    df = pd.DataFrame(valid_rows, columns=headers)
    df["step"] = _normalize_step_numbers(
        df["step"].astype(int).tolist(), "energyprofile.txt"
    )

    return df


# In[117]:


def plot_pepad(
    df: pd.DataFrame,
    parameters: dict[str, object],
    plot: str | list[str],
    rolling: int,
) -> None:
    """
    Function:
        Plot selected PepAD energy-profile terms against MC step.

    Parameters
    ----------
    df : pandas.DataFrame
        Energy-profile data returned by read_energy_profile().

    parameters : dict
        PepAD parameters including BASE_DIR and N_STEPS.

    plot : str or list of str
        Terms to plot: score, gbind, pagg, or rmsd.

    rolling : int
        Rolling-average window. Use 0 or 1 for unaveraged values.

    Returns
    -------
    None
        This function does not return a value.

    Outputs
    -------
    step_evolution.png
        Selected profile curves written at 600 dpi and displayed.
    """
    # Validate the plotting options.
    if rolling < 0:
        raise ValueError("--rolling must be 0 or a positive integer.")

    if rolling > len(df):
        raise ValueError(
            f"--rolling cannot exceed the {len(df)} rows in energyprofile.txt."
        )

    if isinstance(plot, str):
        plot = [plot]

    plot_items = [item.lower() for item in plot]
    valid_items = ["score", "rmsd", "gbind", "pagg"]
    invalid_items = [item for item in plot_items if item not in valid_items]

    if invalid_items:
        raise ValueError(f"Unknown plot option(s): {', '.join(invalid_items)}")

    colors = {
        "score": "#DB594B",
        "gbind": "#7E57C2",
        "pagg": "#3FA34D",
        "rmsd": "#4B89DB",
    }

    plot_columns = {
        "score": ("Score", "Score (kcal/mol)"),
        "gbind": ("G_bind", "G_bind (kcal/mol)"),
        "pagg": ("Pagg", "Pagg (kcal/mol)"),
        "rmsd": ("rmsd", "RMSD (\u212B)"),
    }

    # Plot the selected terms against the MC step.
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax2 = ax.twinx() if "rmsd" in plot_items else None

    for plot_name in plot_items:
        column, label = plot_columns[plot_name]

        if rolling > 1:
            y = df[column].rolling(rolling).mean()
        else:
            y = df[column]

        plot_axis = ax2 if plot_name == "rmsd" else ax
        plot_axis.plot(
            df["step"],
            y,
            label=label,
            color=colors[plot_name],
            linewidth=1,
        )

    ax.set_xlabel("Step", fontsize=18)
    ax.set_ylabel("Energy (kcal/mol)", fontsize=20)
    ax.tick_params(axis="both", labelsize=18)
    final_step = int(df["step"].max()) if not df.empty else 0
    ax.set_xlim(0, max(int(parameters["N_STEPS"]), final_step))

    lines, labels = ax.get_legend_handles_labels()

    if ax2 is not None:
        ax2.set_ylabel("RMSD (\u212B)", fontsize=20)
        ax2.tick_params(axis="both", labelsize=18)
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        lines = lines + lines_2
        labels = labels + labels_2

    ax.legend(
        lines,
        labels,
        loc="lower center",
        ncol=len(plot_items),
        bbox_to_anchor=(0.5, 1.01),
        fontsize=16,
        frameon=False,
    )

    fig.savefig(
        os.path.join(parameters["BASE_DIR"], "step_evolution.png"),
        dpi=600,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(fig)
    print("Saved: step_evolution.png")


# In[118]:


def generate_pepad_report(
    df: pd.DataFrame, parameters: dict[str, object], top: int = 10
) -> pd.DataFrame:
    """
    Function:
        Select the best unique sequences and write PepAD report.txt.

    Parameters
    ----------
    df : pandas.DataFrame
        Energy-profile data returned by read_energy_profile().

    parameters : dict
        PepAD parameters including BASE_DIR.

    top : int, default=10
        Number of unique sequences.

    Returns
    -------
    pandas.DataFrame
        Best unique sequences in ascending Score order.

    Outputs
    -------
    PepAD report.txt
        Table of the selected unique peptide sequences.
    """
    # Validate the report option.
    if top < 1:
        raise ValueError("--top must be greater than 0.")

    # Rank the best unique peptide sequences.
    sequence_counts = df["Sequence"].value_counts().rename("Counts")
    df_sorted = df.sort_values(by="Score", ascending=True)
    df_unique_best = df_sorted.drop_duplicates(subset="Sequence", keep="first")
    df_unique_best = df_unique_best.merge(
        sequence_counts, how="left", left_on="Sequence", right_index=True
    )
    df_min_unique = df_unique_best.head(top)

    # Write the PepAD report.
    filename = os.path.join(parameters["BASE_DIR"], "PepAD report.txt")
    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"---{top} unique peptides with best score (energy profile)---\n")
        f.write(df_min_unique.to_string(index=False, float_format="%.2f"))
        f.write("\n\n")

    # print("Saved: PepAD report.txt")
    # print(df_min_unique.to_string(index=False, float_format="%.2f"))

    return df_min_unique


# In[ ]:


def get_parser() -> argparse.ArgumentParser:
    """
    Function:
        Define terminal flags for profile and energy-detail analysis.

    Parameters
    ----------
    None

    Returns
    -------
    argparse.ArgumentParser
        Configured PepAD analyzer command-line parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "PepAD analyzer reads parameters from the initial-structure PDB "
            "file and input.txt, then analyzes energyprofile.txt and energydetails.txt."
        )
    )
    parser.add_argument(
        "--directory",
        default=".",
        help="PepAD run directory; default is the current directory",
    )
    parser.add_argument(
        "--profiles",
        action="store_true",
        help="Analyze energyprofile.txt and output unique peptides to PepAD_report",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Integer. Number of unique peptides written to PepAD_report",
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="Analyze energydetails.txt and write Detail_report",
    )
    parser.add_argument(
        "--plot",
        nargs="+",
        choices=["score", "rmsd", "gbind", "pagg"],
        default=[],
        help="Optional step-wise plots: [score] [rmsd] [gbind] [pagg]",
    )
    parser.add_argument(
        "--violin",
        nargs="*",
        choices=["energy", "contribution"],
        default=None,
        help=(
            "Optional detail violin plots: [energy] [contribution]. "
            "Without a value, energy is plotted."
        ),
    )
    parser.add_argument(
        "--rolling",
        type=int,
        default=0,
        help="Integer. Number of rolling-average windows for plots",
    )
    return parser


# In[120]:


def read_energy_detail(
    basepath: str, expected_sequence_length: int | None = None
) -> pd.DataFrame:
    """
    Function:
        Read energy details, including trials rejected by the rotamer check.

    Parameters
    ----------
    basepath : str
        Directory containing energydetails.txt.

    expected_sequence_length : int or None, default=None
        Expected residue count including caps. The first record defines it
        when omitted.

    Returns
    -------
    pandas.DataFrame
        Trial energies, delta energies, outcome, and peptide sequence.
    """
    energy_file = os.path.join(basepath, "energydetails.txt")
    headers = [
        "step",
        "attempt",
        "Score",
        "E_VDW",
        "E_ELE",
        "E_SGB",
        "E_SUR",
        "TS",
        "PAGG",
        "dE_VDW",
        "dE_ELE",
        "dE_SGB",
        "dE_SUR",
        "dTS",
        "dPAGG",
        "d_abs_Score",
        "type",
        "acceptance",
        "Sequence",
    ]

    if not os.path.isfile(energy_file):
        raise FileNotFoundError(f"energydetails.txt is not found: {energy_file}")

    with open(energy_file, "r", encoding="utf-8", errors="replace") as f:
        data_lines = [
            (line_number, line.split())
            for line_number, line in enumerate(f, start=1)
            if line.strip()
        ]

    if not data_lines or len(data_lines) % 2 != 0:
        raise ValueError(
            "Energy detail is incomplete. Each record must contain an energy "
            "line and a sequence line."
        )

    rows = []
    sequence_length = expected_sequence_length
    valid_outcomes = {"Accept", "Reject-MC", "Reject-Rotamer"}
    for i in range(0, len(data_lines), 2):
        line_number, values = data_lines[i]
        sequence_line, sequence = data_lines[i + 1]

        if len(values) not in (4, 11, 17, 18):
            raise ValueError(
                "Energy detail is incomplete. Missing or extra value at line "
                f"{line_number}."
            )
        if values[-1] not in valid_outcomes:
            raise ValueError(
                "Energy detail is incomplete. Missing or invalid trial "
                f"outcome at line {line_number}."
            )

        sequence_values = [value.upper() for value in sequence]
        if sequence_length is None:
            sequence_length = len(sequence_values)
        sequence_values = _validate_residue_sequence(
            sequence_values,
            sequence_length,
            f"energydetails.txt line {sequence_line}",
        )

        try:
            if len(values) == 4:
                if values[3] != "Reject-Rotamer":
                    raise ValueError
                energy_values = [int(values[0]), int(values[1])] + [None] * 7
                delta_values = [None] * 7
            else:
                energy_values = [int(values[0]), int(values[1])] + [
                    float(value) for value in values[2:9]
                ]
                if len(values) == 17:
                    delta_values = [float(value) for value in values[9:15]] + [None]
                elif len(values) == 18:
                    delta_values = [float(value) for value in values[9:16]]
                else:
                    delta_values = [None] * 7
        except ValueError:
            raise ValueError(
                f"Energy detail is incomplete. Invalid number at line {line_number}."
            )

        rows.append(
            energy_values
            + delta_values
            + [values[-2], values[-1], "-".join(sequence_values)]
        )

    energy_detail = pd.DataFrame(rows, columns=headers)
    energy_detail["step"] = _normalize_step_numbers(
        energy_detail["step"].astype(int).tolist(),
        "energydetails.txt",
        energy_detail["attempt"].astype(int).tolist(),
    )

    energy_detail["dE_(ELE+SGB)"] = energy_detail["dE_ELE"] + energy_detail["dE_SGB"]
    energy_detail["d_abs_Score"] = (
        abs(energy_detail["dE_VDW"])
        + abs(energy_detail["dE_(ELE+SGB)"])
        + abs(energy_detail["dE_SUR"])
        + abs(energy_detail["dTS"])
        + abs(energy_detail["dPAGG"])
    )

    energy_detail["dE_VDW_p"] = (
        energy_detail["dE_VDW"] / energy_detail["d_abs_Score"]
    )  # percentage of VDW changes
    energy_detail["dE_(ELE+SGB)_p"] = (
        energy_detail["dE_(ELE+SGB)"] / energy_detail["d_abs_Score"]
    )  # percentage of ELE changes
    energy_detail["dE_SUR_p"] = (
        energy_detail["dE_SUR"] / energy_detail["d_abs_Score"]
    )  # percentage of non-polar solvation changes
    energy_detail["dTS_p"] = (
        energy_detail["dTS"] / energy_detail["d_abs_Score"]
    )  # percentage of TS changes
    energy_detail["dPAGG_p"] = (
        energy_detail["dPAGG"] / energy_detail["d_abs_Score"]
    )  # percentage of TS changes

    trial_counts = energy_detail["acceptance"].value_counts()
    print(f"Read {len(energy_detail)} energy-detail records.")
    print("Trial counts:")
    print(f"  Accept: {trial_counts.get('Accept', 0)}")
    print(f"  Reject-MC: {trial_counts.get('Reject-MC', 0)}")
    print(f"  Reject-Rotamer: {trial_counts.get('Reject-Rotamer', 0)}")
    return energy_detail


# In[121]:


def plot_energy_detail(energy_detail: pd.DataFrame) -> None:
    """
    Function:
        Plot energy values and their trial changes against MC step.

    Parameters
    ----------
    energy_detail : pandas.DataFrame
        Energy_detail data returned by read_energy_detail().

    Returns
    -------
    None
        This function does not return a value.

    Outputs
    -------
    matplotlib figure
        Energy and delta-energy curves displayed in the notebook.
    """
    terms = [
        ("Score", "d_abs_Score"),
        ("E_VDW", "dE_VDW"),
        ("E_(ELE+SGB)", "dE_(ELE+SGB)"),
        ("E_SUR", "dE_SUR"),
        ("TS", "dTS"),
        ("PAGG", "dPAGG"),
    ]
    energy_color = "#4B89DB"
    delta_color = "#DB594B"
    has_delta = energy_detail["dE_VDW"].notna().any()
    fig, axes = plt.subplots(4, 2, figsize=(16, 18), constrained_layout=True)

    for ax, (energy, delta) in zip(axes.flat, terms):
        ax.plot(
            energy_detail["step"],
            energy_detail[energy],
            color=energy_color,
            linewidth=1,
        )
        ax.set_title(energy, fontsize=18)
        ax.set_xlabel("Step", fontsize=18)
        ax.set_ylabel(energy, fontsize=18, color=energy_color)
        ax.tick_params(axis="both", labelsize=16)
        ax.tick_params(axis="y", colors=energy_color)

        delta_axis = ax.twinx()
        if energy_detail[delta].notna().any():
            delta_axis.plot(
                energy_detail["step"],
                energy_detail[delta],
                color=delta_color,
                linewidth=1,
            )
            delta_axis.set_ylabel(delta, fontsize=18, color=delta_color)
            delta_axis.tick_params(axis="y", labelsize=16, colors=delta_color)
        else:
            delta_axis.set_visible(False)

    fig.delaxes(axes[3, 1])
    if not has_delta:
        print("This energydetails.txt file does not contain delta-energy columns.")
    plt.show()


# In[122]:


# The complete notebook example is after analyze_delta_energy().


# In[ ]:


def analyze_delta_energy(
    energy_detail: pd.DataFrame,
    output_file: str = "Detail_report.txt",
    plot_options: str | list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Function:
        Analyze how each energy term contributes to accepted and rejected
        trial moves, draw percentage and raw-value violin plots, and write
        the detailed trial report. Raw plots use negative accepted values
        and positive Reject-MC values.

    Parameters
    ----------
    energy_detail : pandas.DataFrame
        Trial data returned by read_energy_detail().

    output_file : str, default='Detail_report.txt'
        Path for the energy-detail report.

    plot_options : list of str or None
        Optional detail plots: energy, contribution, or both.

    Equations
    ---------
    Delta_Ei is one of Delta_VDW, Delta_ELE_SGB, Delta_SUR,
    Delta_(-TS), or Delta_(-PAGG).
    Delta_ELE_SGB = Delta_ELE + Delta_SGB
    Delta_(-TS) = -Delta_TS
    Delta_(-PAGG) = -Delta_PAGG
    Delta_Score = sum(Delta_Ei)
    Total_negative = sum(abs(min(Delta_Ei, 0)))
    Total_positive = sum(max(Delta_Ei, 0))
    Negative_Ei_percent = abs(min(Delta_Ei, 0)) / Total_negative * 100
    Positive_Ei_percent = max(Delta_Ei, 0) / Total_positive * 100
    PepAD accepts Delta_Score <= 0 directly. For Delta_Score > 0,
    acceptance probability = exp(-Delta_Score / kBT_seq).

    Returns
    -------
    tuple of pandas.DataFrame
        Trial-level dominant terms, dominant counts, and percentages.

    Outputs
    -------
    Detail_report.txt
        Energy-detail report written whenever this function runs.

    Delta_energy_distribution.png
        Raw energy violin plot written at 600 dpi when energy is selected.

    Delta_energy_contribution.png
        Percentage violin plot written at 600 dpi when contribution is selected.
    """
    terms = ["VDW", "ELE+SGB", "SUR", "(-TS)", "(-Pagg)"]
    colors = {"Accept": "#4B89DB", "Reject-MC": "#DB594B"}
    subplot_pad = 10
    if plot_options is None:
        plot_options = []
    elif isinstance(plot_options, str):
        plot_options = [plot_options]
    plot_options = [option.lower() for option in plot_options]
    invalid_plots = set(plot_options) - {"energy", "contribution"}
    if invalid_plots:
        invalid_text = ", ".join(sorted(invalid_plots))
        raise ValueError(f"Unknown detail plot option(s): {invalid_text}")
    percentage_ylabel = (
        r"$\mathrm{Accept}:\ "
        r"|\Delta\Delta E_k^{(-)}|\,/\,"
        r"\sum_i |\Delta\Delta E_i^{(-)}|\times 100\%$"
        "\n"
        r"$\mathrm{Reject\!\!-\!\!MC}:\ "
        r"|\Delta\Delta E_k^{(+)}|\,/\,"
        r"\sum_i |\Delta\Delta E_i^{(+)}|\times 100\%$"
    )
    raw_ylabel = (
        r"$\mathrm{Accept}:\ \Delta\Delta E_k^{(-)}\ "
        r"(\mathrm{kcal/mol})$"
        "\n"
        r"$\mathrm{Reject\!\!-\!\!MC}:\ "
        r"\Delta\Delta E_k^{(+)}\ (\mathrm{kcal/mol})$"
    )

    # These are signed contributions to Score. TS and PAGG enter Score with minus signs.
    contributions = pd.DataFrame(
        {
            "VDW": energy_detail["dE_VDW"],
            "ELE+SGB": energy_detail["dE_(ELE+SGB)"],
            # "SGB": energy_detail["dE_SGB"],
            "SUR": energy_detail["dE_SUR"],
            "(-TS)": -energy_detail["dTS"],
            "(-Pagg)": -energy_detail["dPAGG"],
        }
    )
    complete = contributions.notna().all(axis=1)
    accepted = complete & energy_detail["acceptance"].eq("Accept")
    rejected = complete & energy_detail["acceptance"].eq("Reject-MC")
    if not accepted.any() and not rejected.any():
        raise ValueError(
            "No delta-energy data are available for Accept or Reject-MC trials."
        )

    negative = (
        contributions.loc[accepted].clip(upper=0).abs()
    )  # pickout negative values in VMD ELE+SGB SUR (-TS) (-PAGG), and assign as accept
    positive = contributions.loc[rejected].clip(
        lower=0
    )  # pickout positive values in VMD ELE+SGB SUR (-TS) (-PAGG), and assign as reject
    negative_total = negative.sum(axis=1)  # negative contribution to score change
    positive_total = positive.sum(axis=1)  # positive contribution to score change
    negative_percentage = (
        negative.loc[negative_total > 0]
        .div(negative_total[negative_total > 0], axis=0)
        .mul(100)
    )
    positive_percentage = (
        positive.loc[positive_total > 0]
        .div(positive_total[positive_total > 0], axis=0)
        .mul(100)
    )
    accepted_signed = contributions.loc[accepted]
    rejected_signed = contributions.loc[rejected]
    accepted_negative_median = accepted_signed.where(
        accepted_signed < 0
    ).median()
    rejected_positive_median = rejected_signed.where(
        rejected_signed > 0
    ).median()

    contribution_summary = pd.DataFrame(
        {
            "Median_negative_ddEi_when_accepted(kcal/mol)": (
                accepted_negative_median
            ),
            "Median_positive_ddEi_when_rejected(kcal/mol)": (
                rejected_positive_median
            ),
            "Frequency_of_negative_ddEi_when_accepted(%)": (
                negative.gt(0).mean().mul(100)
            ),
            "Median_negative_ddEi_contribution_when_accepted(%)": (
                negative_percentage.mask(negative_percentage.eq(0)).median()
            ),
            "Frequency_of_positive_ddEi_when_rejected(%)": (
                positive.gt(0).mean().mul(100)
            ),
            "Median_positive_ddEi_contribution_when_rejected(%)": (
                positive_percentage.mask(positive_percentage.eq(0)).median()
            ),
        }
    ).round(2)

    output_directory = os.path.dirname(os.path.abspath(output_file))

    ################################
    # Plot Delta_Ei distribution
    ################################
    if "energy" in plot_options:
        fig_raw, axes_raw = plt.subplots(
            2, 3, figsize=(18, 11), sharey=False, constrained_layout=True
        )
        fig_raw.set_constrained_layout_pads(
            w_pad=subplot_pad / 72,
            h_pad=subplot_pad / 72,
        )
        for ax, energy in zip(axes_raw.flat, terms):
            accepted_values = contributions.loc[accepted, energy]
            accepted_values = accepted_values[accepted_values < 0]
            rejected_values = contributions.loc[rejected, energy]
            rejected_values = rejected_values[rejected_values > 0]
            groups = [
                (1, "Accept", accepted_values),
                (2, "Reject-MC", rejected_values),
            ]
            for position, outcome, values in groups:
                values = values.dropna()
                if len(values) > 1 and values.nunique() > 1:
                    violin = ax.violinplot(
                        values,
                        positions=[position],
                        widths=0.75,
                        showmeans=False,
                        showmedians=False,
                        showextrema=False,
                    )
                    violin["bodies"][0].set_facecolor(colors[outcome])
                    violin["bodies"][0].set_edgecolor(colors[outcome])
                    violin["bodies"][0].set_alpha(0.65)
                elif len(values) > 0:
                    ax.scatter(
                        position, values.median(), color=colors[outcome], s=45
                    )

                if len(values) > 0:
                    median = values.median()
                    ax.hlines(
                        median,
                        position - 0.18,
                        position + 0.18,
                        color="black",
                        linewidth=2.5,
                        zorder=4,
                    )

            ax.axhline(0, color="gray", linewidth=1, linestyle="--")
            ax.set_title(energy, fontsize=18, pad=10)
            ax.set_xlim(0.5, 2.5)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Accept\nnegative", "Reject-MC\npositive"])
            ax.tick_params(axis="both", labelsize=16)

        axes_raw.flat[-1].axis("off")
        fig_raw.supylabel(raw_ylabel, fontsize=17)
        energy_figure = os.path.join(
            output_directory, "Delta_energy_distribution.png"
        )
        fig_raw.savefig(energy_figure, dpi=600, bbox_inches="tight")
        plt.show()
        plt.close(fig_raw)
        print(f"Energy violin written to: {energy_figure}")

    ################################
    # Plot Delta_Ei percentages
    ################################
    if "contribution" in plot_options:
        fig, axes = plt.subplots(
            2, 3, figsize=(18, 11), sharey=True, constrained_layout=True
        )
        fig.set_constrained_layout_pads(
            w_pad=subplot_pad / 72,
            h_pad=subplot_pad / 72,
        )
        for ax, energy in zip(axes.flat, terms):
            groups = [
                (1, "Accept", negative_percentage[energy]),
                (2, "Reject-MC", positive_percentage[energy]),
            ]
            for position, outcome, values in groups:
                values = values[values > 0].dropna()
                if len(values) > 1 and values.nunique() > 1:
                    violin = ax.violinplot(
                        values,
                        positions=[position],
                        widths=0.75,
                        showmeans=False,
                        showmedians=False,
                        showextrema=False,
                    )
                    violin["bodies"][0].set_facecolor(colors[outcome])
                    violin["bodies"][0].set_edgecolor(colors[outcome])
                    violin["bodies"][0].set_alpha(0.65)
                elif len(values) > 0:
                    ax.scatter(
                        position, values.median(), color=colors[outcome], s=45
                    )

                if len(values) > 0:
                    median = values.median()
                    ax.hlines(
                        median,
                        position - 0.18,
                        position + 0.18,
                        color="black",
                        linewidth=2.5,
                        zorder=4,
                    )
                    ax.text(
                        position + 0.20,
                        median,
                        f"{median:.1f}",
                        fontsize=12,
                        va="center",
                    )

            ax.set_title(energy, fontsize=18, pad=10)
            ax.set_xlim(0.5, 2.5)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(
                [
                    f"Accept\n{negative[energy].gt(0).mean() * 100:.1f}% negative",
                    f"Reject-MC\n{positive[energy].gt(0).mean() * 100:.1f}% positive",
                ]
            )
            ax.set_ylim(0, 100)
            ax.tick_params(axis="both", labelsize=16)

        axes.flat[-1].axis("off")
        fig.supylabel(percentage_ylabel, fontsize=17)
        percent_figure = os.path.join(
            output_directory, "Delta_energy_contribution.png"
        )
        fig.savefig(percent_figure, dpi=600, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Contribution violin written to: {percent_figure}")

    accepted_trials = energy_detail.loc[
        negative_percentage.index, ["step", "attempt", "type", "acceptance", "Sequence"]
    ].copy()
    accepted_trials["group"] = "Accept: negative"
    accepted_trials["dominant_energy"] = negative_percentage.idxmax(axis=1)
    accepted_trials["highest_percentage"] = negative_percentage.max(axis=1).round(2)
    accepted_trials["selected_total"] = negative_total.loc[negative_percentage.index]
    accepted_trials["delta_score"] = contributions.loc[negative_percentage.index].sum(
        axis=1
    )

    rejected_trials = energy_detail.loc[
        positive_percentage.index, ["step", "attempt", "type", "acceptance", "Sequence"]
    ].copy()
    rejected_trials["group"] = "Reject-MC: positive"
    rejected_trials["dominant_energy"] = positive_percentage.idxmax(axis=1)
    rejected_trials["highest_percentage"] = positive_percentage.max(axis=1).round(2)
    rejected_trials["selected_total"] = positive_total.loc[positive_percentage.index]
    rejected_trials["delta_score"] = contributions.loc[positive_percentage.index].sum(
        axis=1
    )

    dominant_trials = pd.concat([accepted_trials, rejected_trials]).sort_index()
    dominant_counts = pd.crosstab(
        dominant_trials["dominant_energy"], dominant_trials["group"]
    )
    dominant_counts = dominant_counts.reindex(
        index=terms, columns=["Accept: negative", "Reject-MC: positive"], fill_value=0
    )
    dominant_percentage = (
        dominant_counts.div(dominant_counts.sum(axis=0), axis=1).mul(100).round(1)
    )

    dominant_summary = pd.DataFrame(
        {
            "Dominant_energy": terms,
            "Accepts": dominant_counts["Accept: negative"].to_numpy(),
            "Accept_rate(%)": dominant_percentage["Accept: negative"].to_numpy(),
            "Rejects-MC": dominant_counts["Reject-MC: positive"].to_numpy(),
            "Rejects-MC_rate(%)": dominant_percentage["Reject-MC: positive"].to_numpy(),
        }
    )
    dominant_summary.loc[len(dominant_summary)] = [
        "Total",
        int(dominant_counts["Accept: negative"].sum()),
        100.0,
        int(dominant_counts["Reject-MC: positive"].sum()),
        100.0,
    ]
    trial_counts = energy_detail["acceptance"].value_counts()
    accept_count = int(trial_counts.get("Accept", 0))
    reject_mc_count = int(trial_counts.get("Reject-MC", 0))
    reject_rotamer_count = int(trial_counts.get("Reject-Rotamer", 0))
    accepted_without_negative = int((negative_total == 0).sum())
    rejected_without_positive = int((positive_total == 0).sum())
    accepted_uphill = int((contributions.loc[accepted].sum(axis=1) > 0).sum())

    print(
        f"Accepted trials without a negative contribution: {accepted_without_negative}"
    )
    print(
        f"Reject-MC trials without a positive contribution: {rejected_without_positive}"
    )
    print(f"Accepted uphill trials (delta Score > 0): {accepted_uphill}")

    #############################
    # Output Detail_report.txt
    #############################
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f"Total number of trials: {len(energy_detail)}\n")
        f.write(f"Accepted trials: {accept_count}\n")
        f.write(f"Reject-MC trials: {reject_mc_count}\n")
        f.write(f"Reject-Rotamer trials: {reject_rotamer_count}\n")
        f.write(
            "Accepted trials without a negative contribution: "
            f"{accepted_without_negative}\n"
        )
        f.write(
            "Reject-MC trials without a positive contribution: "
            f"{rejected_without_positive}\n"
        )
        f.write(f"Accepted uphill trials (delta Score > 0): {accepted_uphill}\n\n")
        f.write(
            "Energy and percentage contributions to accepted and rejected "
            "trials:\n"
        )
        f.write(
            "Energy_term\tMedian_negative_ddEi_when_accepted(kcal/mol)\t"
            "Median_positive_ddEi_when_rejected(kcal/mol)\t"
            "Frequency_of_negative_ddEi_when_accepted(%)\t"
            "Median_negative_ddEi_contribution_when_accepted(%)\t"
            "Frequency_of_positive_ddEi_when_rejected(%)\t"
            "Median_positive_ddEi_contribution_when_rejected(%)\n"
        )
        for energy, row in contribution_summary.iterrows():
            values = [
                "NA" if pd.isna(value) else f"{value:.2f}" for value in row
            ]
            f.write(energy + "\t" + "\t".join(values) + "\n")
        f.write("\n\n")
        f.write("Dominant contribution count and percentage:\n")
        f.write(
            "Dominant_energy\tAccepts\tAccept(%)\t"
            "Rejects-MC\tRejects-MC(%)\n"
        )
        for row in dominant_summary.itertuples(index=False, name=None):
            f.write(
                f"{row[0]}\t{int(row[1])}\t{row[2]:.1f}\t"
                f"{int(row[3])}\t{row[4]:.1f}\n"
            )

    print(f"Detail report written to: {output_file}")
    return dominant_trials, dominant_counts, dominant_percentage


# In[ ]:


# Notebook example: run every profile and energy-detail analysis.
# Uncomment this block to run the complete example.
# parameters = read_input(".")
# base_dir = parameters["BASE_DIR"]

# # Plot every energy-profile quantity and write the top-10 report.
# profile_length = parameters["N_AA"]
# energy_profile = read_energy_profile(base_dir, profile_length)
# profile_plots = ["score", "rmsd", "gbind", "pagg"]
# plot_pepad(energy_profile, parameters, profile_plots, rolling=0)
# best_peptides = generate_pepad_report(
#     energy_profile, parameters, top=10
# )

# # Plot every energy-detail trace and both violin plot types.
# detail_length = parameters["N_RESIDUES"]
# energy_detail = read_energy_detail(base_dir, detail_length)
# plot_energy_detail(energy_detail)
# detail_report = os.path.join(base_dir, "Detail_report.txt")
# detail_results = analyze_delta_energy(
#     energy_detail, detail_report, ["energy", "contribution"]
# )
# dominant_trials, dominant_counts, dominant_percentage = detail_results


# In[ ]:


def main(argv: list[str] | None = None) -> None:
    """
    Function:
        Run profile analysis, energy-detail analysis, or both.

    Parameters
    ----------
    argv : list of str or None, default=None
        Terminal arguments. None reads the real command line.

    Returns
    -------
    None
        This function does not return a value.

    Outputs
    -------
    PepAD report.txt
        Written when --profiles is selected.

    Detail_report.txt
        Written when --details is selected.

    PNG figures
        Written only when --plot or --violin is selected.
    """
    parser = get_parser()
    args = parser.parse_args(argv)
    if not args.profiles and not args.details:
        parser.print_help()
        return

    profile_plots = list(dict.fromkeys(args.plot))
    if args.violin is None:
        detail_plots = []
    elif not args.violin:
        detail_plots = ["energy"]
    else:
        detail_plots = [
            option
            for option in ("energy", "contribution")
            if option in args.violin
        ]

    if profile_plots and not args.profiles:
        parser.error("Profile plot options require --profiles.")
    if detail_plots and not args.details:
        parser.error("Detail violin plot options require --details.")

    parameters = read_input(args.directory)
    if parameters is None:
        return

    if args.profiles:
        profile = read_energy_profile(
            parameters["BASE_DIR"], parameters["N_AA"]
        )
        if profile_plots:
            plot_pepad(profile, parameters, profile_plots, args.rolling)
        generate_pepad_report(profile, parameters, args.top)

    if args.details:
        energy_detail = read_energy_detail(
            parameters["BASE_DIR"], parameters["N_RESIDUES"]
        )
        detail_report = os.path.join(
            parameters["BASE_DIR"], "Detail_report.txt"
        )
        analyze_delta_energy(energy_detail, detail_report, detail_plots)


def cli() -> None:
    """Run the installed ``analyzer`` terminal command."""
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"Analyzer error: {error}") from error


if __name__ == "__main__":
    cli()
