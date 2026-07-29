# Initial Structure Builder v1.23

# Source version: v1.23


# %% Cell 30ee6987
import numpy as np
import pandas as pd
import argparse
import random
from typing import List, Optional, Sequence, Tuple, Union
import math, warnings
from pathlib import Path


DEFAULT_ROTAMER_DIR = Path(__file__).resolve().parent / "RotamerLibrary"


# %% Cell 79accf39-9c01-4f68-ae10-575e04ccb1bd
def linear_fitting(data):
    #      reference:
    # https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d

    """Function:
        Fit one straight line through a group of 3D points and obtain the main line
        direction.

    Parameters
    ----------
    data : numpy.ndarray or pandas.DataFrame with shape (n, 3)
        Each row is one 3D point.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Fitted direction vector [x, y, z].
    """
    datamean = data.mean(axis=0)
    uu, dd, vv = np.linalg.svd(data - datamean)
    return vv[0]


def unit_vector(vector, label="vector", allow_zero=False):
    """Function:
        Normalize a vector so its length becomes 1.

    Parameters
    ----------
    vector : list, tuple, or numpy.ndarray
        Vector components to normalize.
    label : str
        Name included in an error message when the vector is invalid.
    allow_zero : bool
        When True, a near-zero vector is allowed and returns None.

    Returns
    -------
    numpy.ndarray or None
        Array with the same shape as vector, or None.
    """
    vector = np.asarray(vector, dtype=np.float64)
    norm = np.linalg.norm(vector)
    if norm < 1e-8:
        if allow_zero:
            return None
        raise ValueError(f"Cannot normalize near-zero {label}.")
    return vector / norm


def cross2d(a, b):
    """Function:
        Calculate the signed 2D cross product used to judge clockwise or
        counterclockwise rotation.

    Parameters
    ----------
    a : array-like
        First two-dimensional vector.
    b : array-like
        Second two-dimensional vector.

    Returns
    -------
    float
        Signed scalar cross product a_x*b_y - a_y*b_x.
    """
    return a[0] * b[1] - a[1] * b[0]


def rotation_angle(v1, v2):
    """Function:
        Calculate the signed angle needed to rotate one 2D vector onto another vector.

    Parameters
    ----------
    v1 : Vector, Atom, or array-like
        Starting two-dimensional vector.
    v2 : Vector, Atom, or array-like
        Target two-dimensional vector.

    Returns
    -------
    float
        Signed rotation angle in radians.
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    cos_angle = np.dot(v1_u, v2_u) / (
        np.sqrt(np.dot(v1_u, v1_u)) * np.sqrt(np.dot(v2_u, v2_u))
    )
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    angle = np.arccos(cos_angle)

    # judge rotation direction
    cross_product = cross2d(
        v1_u, v2_u
    )  # cross_product OF 2d vector, it return a scalar
    if cross_product < 0:
        direction = -1
    else:
        direction = 1

    angle = angle * direction
    return angle


def rotation_coordinates(df, axis, angle):
    """Function:
        Rotate all atom coordinates in a peptide table around the x, y, or z axis.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    axis : str
        'x', 'y', or 'z' selects the rotation axis and 'none' returns an unchanged copy.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    pandas.DataFrame
        Rotated copy of df.
    """
    df_temp = df.copy()
    data = df_temp[["x", "y", "z"]].to_numpy()

    if axis == "x":
        M = np.array(
            [
                [1, 0, 0],
                [0, np.cos(angle), -np.sin(angle)],
                [0, np.sin(angle), np.cos(angle)],
            ]
        )
    elif axis == "y":
        M = np.array(
            [
                [np.cos(angle), 0, np.sin(angle)],
                [0, 1, 0],
                [-np.sin(angle), 0, np.cos(angle)],
            ]
        )
    elif axis == "z":
        M = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )
    elif axis == "none":
        return df_temp

    rotated_data = (M @ data.T).T  # Transpose to apply matrix multiplication correctly
    df_temp[["x", "y", "z"]] = rotated_data

    return df_temp


def add_NH_remove_OH(peptide_df: pd.DataFrame) -> pd.DataFrame:
    """Function:
        Prepare uncapped peptide termini directly in a peptide DataFrame.

    Parameters
    ----------
    peptide_df : pandas.DataFrame
        Generated peptide atoms before terminal processing.

    Returns
    -------
    pandas.DataFrame
        Peptide atoms with terminal hydrogens and residue names corrected.
    """
    df = peptide_df.copy().reset_index(drop=True)
    ##### remove "H" on "OXT" ###############
    df = df[df["atom_name"] != "HXT"]
    df.loc[(df["atom_name"] == "H") & (df["resid"] == 1), "atom_name"] = "H1"

    h2 = df[df["atom_name"] == "H2"].iloc[0].copy()
    h3 = h2.copy()
    h3["atom_name"] = "H3"
    idx = df[df["atom_name"] == "H2"].index[0] + 1  # index of H2 in the df
    df = pd.concat([df.iloc[:idx], pd.DataFrame([h3]), df.iloc[idx:]]).reset_index(
        drop=True
    )  # combine dataframe
    df["anum"] = range(1, len(df) + 1)  # Reindex atom number

    ################## calculate H coordinates #########
    # Solve the corrected system of equations
    if df[df["resid"] == 1].iloc[0]["aa_name"] == "GLY":
        cb = (
            df[(df["atom_name"] == "HA2") & (df["resid"] == 1)]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
    else:
        cb = (
            df[(df["atom_name"] == "CB") & (df["resid"] == 1)]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )

    ca = (
        df[(df["atom_name"] == "CA") & (df["resid"] == 1)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )
    n = (
        df[(df["atom_name"] == "N") & (df["resid"] == 1)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )

    # need to reconstruct the numpy array again to prevent error
    cb = np.array([cb[0], cb[1], cb[2]])
    ca = np.array([ca[0], ca[1], ca[2]])
    n = np.array([n[0], n[1], n[2]])

    # Parameters
    dihedral_angle_rad = 59.98 / 180 * np.pi
    bond_angle_rad = 109.49 / 180 * np.pi
    dihedral_angle = 59.98
    bond_angle = 109.49
    bond_length = 1.01

    h_coordinate = place_fourth_atom(cb, ca, n, bond_length, bond_angle, dihedral_angle)

    angles = [120, 240]  # Angles to rotate by
    rotated_h_coords = rotate_around_axis(ca, n, h_coordinate, angles)

    df.loc[df["atom_name"] == "H1", ["x", "y", "z"]] = h_coordinate
    df.loc[df["atom_name"] == "H2", ["x", "y", "z"]] = rotated_h_coords[0]
    df.loc[df["atom_name"] == "H3", ["x", "y", "z"]] = rotated_h_coords[1]

    first_resid = df["resid"].min()
    last_resid = df["resid"].max()
    first_mask = df["resid"] == first_resid
    last_mask = (df["resid"] == last_resid) & ~first_mask
    df.loc[first_mask, "aa_name"] = (
        "N" + df.loc[first_mask, "aa_name"].astype(str)
    )
    df.loc[last_mask, "aa_name"] = (
        "C" + df.loc[last_mask, "aa_name"].astype(str)
    )
    return df.reset_index(drop=True)

def add_caps(peptide_df: pd.DataFrame, type_flag: int) -> pd.DataFrame:
    """Function:
        Add ACE and the selected C-terminal cap directly to a peptide DataFrame.

    Parameters
    ----------
    peptide_df : pandas.DataFrame
        Generated peptide atoms before terminal processing.
    type_flag : int
        1 adds ACE+NME and 2 adds ACE+NHE.

    Returns
    -------
    pandas.DataFrame
        Peptide atoms with the selected terminal caps.
    """
    NHE_data = [
        {
            "anum": 1,
            "atom_name": "N",
            "aa_name": "NHE",
            "resid": 1,
            "x": 2.194,
            "y": 1.598,
            "z": -0.000,
        },
        {
            "anum": 2,
            "atom_name": "HN1",
            "aa_name": "NHE",
            "resid": 1,
            "x": 3.035,
            "y": 1.039,
            "z": -0.000,
        },
        {
            "anum": 3,
            "atom_name": "HN2",
            "aa_name": "NHE",
            "resid": 1,
            "x": 2.250,
            "y": 2.606,
            "z": -0.000,
        },
    ]

    NME_data = [
        {
            "anum": 1,
            "atom_name": "N",
            "aa_name": "NME",
            "resid": 1,
            "x": 3.326,
            "y": 1.548,
            "z": -0.000,
        },
        {
            "anum": 2,
            "atom_name": "H",
            "aa_name": "NME",
            "resid": 1,
            "x": 3.909,
            "y": 0.724,
            "z": -0.000,
        },
        {
            "anum": 3,
            "atom_name": "CH3",
            "aa_name": "NME",
            "resid": 1,
            "x": 3.970,
            "y": 2.846,
            "z": -0.000,
        },
        {
            "anum": 4,
            "atom_name": "HH31",
            "aa_name": "NME",
            "resid": 1,
            "x": 3.212,
            "y": 3.629,
            "z": 0.000,
        },
        {
            "anum": 5,
            "atom_name": "HH32",
            "aa_name": "NME",
            "resid": 1,
            "x": 4.592,
            "y": 2.943,
            "z": 0.890,
        },
        {
            "anum": 6,
            "atom_name": "HH33",
            "aa_name": "NME",
            "resid": 1,
            "x": 4.592,
            "y": 2.943,
            "z": -0.890,
        },
    ]

    ACE_data = [
        {
            "anum": 1,
            "atom_name": "H1",
            "aa_name": "ACE",
            "resid": 1,
            "x": 2.000,
            "y": 1.000,
            "z": -0.000,
        },
        {
            "anum": 2,
            "atom_name": "CH3",
            "aa_name": "ACE",
            "resid": 1,
            "x": 2.000,
            "y": 2.090,
            "z": 0.000,
        },
        {
            "anum": 3,
            "atom_name": "H2",
            "aa_name": "ACE",
            "resid": 1,
            "x": 1.486,
            "y": 2.454,
            "z": 0.890,
        },
        {
            "anum": 4,
            "atom_name": "H3",
            "aa_name": "ACE",
            "resid": 1,
            "x": 1.486,
            "y": 2.454,
            "z": -0.890,
        },
        {
            "anum": 5,
            "atom_name": "C",
            "aa_name": "ACE",
            "resid": 1,
            "x": 3.427,
            "y": 2.641,
            "z": -0.000,
        },
        {
            "anum": 6,
            "atom_name": "O",
            "aa_name": "ACE",
            "resid": 1,
            "x": 4.391,
            "y": 1.877,
            "z": -0.000,
        },
    ]

    NHE = pd.DataFrame(NHE_data)
    NME = pd.DataFrame(NME_data)
    ACE = pd.DataFrame(ACE_data)

    df = peptide_df.copy().reset_index(drop=True)

    ##### remove "H" on "OXT" ###############

    df = df[(df["atom_name"] != "HXT")]  # remove H on OXT

    # Keep one backbone N-H for the ACE amide bond. Put it at the same
    # geometry as the uncapped N-terminal H3, then orient ACE from that H.
    n_term_h_names = ["H", "H1", "HN", "HN1", "H2", "HN2", "H3", "HN3"]
    first_res = df[df["resid"] == 1]
    n = (
        first_res[first_res["atom_name"] == "N"]
        .iloc[0][["x", "y", "z"]]
        .to_numpy(dtype=np.float64)
    )
    ca = (
        first_res[first_res["atom_name"] == "CA"]
        .iloc[0][["x", "y", "z"]]
        .to_numpy(dtype=np.float64)
    )
    if first_res.iloc[0]["aa_name"] == "GLY":
        ref_atom = first_res[first_res["atom_name"].isin(["HA2", "HA"])]
    else:
        ref_atom = first_res[first_res["atom_name"] == "CB"]
    if ref_atom.empty:
        raise ValueError(
            "Cannot add ACE cap: no N-terminal CB/HA reference atom was "
            "found on residue 1."
        )
    ref = ref_atom.iloc[0][["x", "y", "z"]].to_numpy(dtype=np.float64)
    h_seed = place_fourth_atom(ref, ca, n, 1.01, 109.49, 59.98)
    h3_coord = rotate_around_axis(ca, n, h_seed, [120, 240])[1]

    n_term_h = df[(df["resid"] == 1) & (df["atom_name"].isin(n_term_h_names))]
    if n_term_h.empty:
        keep_h_idx = len(df)
        h_row = first_res[first_res["atom_name"] == "N"].iloc[0].copy()
        h_row["atom_name"] = "H"
        h_row[["x", "y", "z"]] = h3_coord
        df = pd.concat([df, pd.DataFrame([h_row])], ignore_index=True)
    else:
        keep_h_idx = n_term_h.index[0]
        remove_h_idx = n_term_h.index[1:]
        df = df.drop(remove_h_idx)
        df.loc[keep_h_idx, "atom_name"] = "H"
        df.loc[keep_h_idx, ["x", "y", "z"]] = h3_coord

    df["resid"] += 1  # increase residue ID by 1

    ####### add ACE group ###########
    df = pd.concat([ACE, df], ignore_index=True)

    ###### calcualte NME or NHE position ########
    last_resid = df[df["atom_name"] == "C"]["resid"].max()
    C = (
        df[(df["atom_name"] == "C") & (df["resid"] == last_resid)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )
    OXT = (
        df[(df["atom_name"] == "OXT") & (df["resid"] == last_resid)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )
    C_OXT_vector = OXT - C
    C_N_vector = C_OXT_vector * 1.291 / np.linalg.norm(C_OXT_vector)
    N = C + C_N_vector
    df = df[~((df["atom_name"] == "OXT") & (df["resid"] == last_resid))]  # remove OXT

    if type_flag == 1:
        n_pos = NME[NME["atom_name"] == "N"][["x", "y", "z"]].values[0]
        ####### add NME group ###########
        NME["x"] = NME["x"] - n_pos[0] + N[0]
        NME["y"] = NME["y"] - n_pos[1] + N[1]
        NME["z"] = NME["z"] - n_pos[2] + N[2]
        NME["resid"] += last_resid

        df = pd.concat([df, NME], ignore_index=True)
        df["anum"] = range(1, len(df) + 1)

        ####### Rotate NME group #######
        CA = (
            df[(df["atom_name"] == "CA") & (df["resid"] == last_resid)]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        N = (
            df[(df["aa_name"] == "NME") & (df["atom_name"] == "N")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        H = (
            df[(df["aa_name"] == "NME") & (df["atom_name"] == "H")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )  # H from NME
        CH3 = (
            df[(df["aa_name"] == "NME") & (df["atom_name"] == "CH3")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        O = (
            df[(df["resid"] == last_resid) & (df["atom_name"] == "O")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        C = np.asarray(C, dtype=np.float64)
        CA = np.asarray(CA, dtype=np.float64)
        N = np.asarray(N, dtype=np.float64)
        H = np.asarray(H, dtype=np.float64)
        CH3 = np.asarray(CH3, dtype=np.float64)
        O = np.asarray(O, dtype=np.float64)

        v1 = C - N
        v2 = CH3 - N
        v3 = rotate_vector_plane_3D(v1, v2, np.pi * 2 / 3) * np.linalg.norm(
            v2
        )  # rotate v2 to v1 by 120 degree
        v4 = rotate_vector_plane_3D(v1, v2, np.pi * 4 / 3) * np.linalg.norm(
            H - N
        )  # also get position of H

        CH3 = N + v3  # correct position of CH3, N
        H = N + v4

        df.loc[
            (df["atom_name"] == "CH3") & (df["aa_name"] == "NME"), ["x", "y", "z"]
        ] = CH3
        df.loc[(df["atom_name"] == "H") & (df["aa_name"] == "NME"), ["x", "y", "z"]] = H

        ################## calculate H coordinates #########
        p1 = H
        p2 = N
        p3 = CH3

        # Parameters
        dihedral_angle_deg = 97.150
        bond_angle_deg = 110.00
        bond_length = 1.078

        h_coordinate = place_fourth_atom(
            p1, p2, p3, bond_length, bond_angle_deg, dihedral_angle_deg
        )  ## important

        angles = [120, 240]  # Angles to rotate by
        rotated_h_coords = rotate_around_axis(p2, p3, h_coordinate, angles)

        df.loc[
            (df["atom_name"] == "HH31") & (df["aa_name"] == "NME"), ["x", "y", "z"]
        ] = h_coordinate
        df.loc[
            (df["atom_name"] == "HH32") & (df["aa_name"] == "NME"), ["x", "y", "z"]
        ] = rotated_h_coords[0]
        df.loc[
            (df["atom_name"] == "HH33") & (df["aa_name"] == "NME"), ["x", "y", "z"]
        ] = rotated_h_coords[1]

        # #################### rotate C(amino acid) - N(NME) axis to 180 degree
        # ################
        theta = dihedral_angle(O, C, N, CH3)
        rotation_angle = 0 * np.pi - theta

        rotation_axis = N - C
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize axis
        rotation_center = N

        # Copy NME group and rotate its atoms
        NME_atoms = df[df["aa_name"] == "NME"].copy()

        for i, row in NME_atoms.iterrows():
            if row["atom_name"] != "N":
                v = (
                    np.array(row[["x", "y", "z"]].values, dtype=np.float64)
                    - rotation_center
                )  # Translate to origin
                v_rot = rotate_vector(v, rotation_axis, rotation_angle)  # Rotate
                NME_atoms.loc[i, ["x", "y", "z"]] = (
                    v_rot + rotation_center
                ).tolist()  # Translate back

        # Replace old NME with rotated version
        df = df[df["aa_name"] != "NME"]
        df = pd.concat([df, NME_atoms], ignore_index=True)

        # Keep terminal O fixed; the C-terminal cap N replaces the removed OXT.

    elif type_flag == 2:
        n_pos = NHE[NHE["atom_name"] == "N"][["x", "y", "z"]].values[0]
        ####### add NHE group ###########
        NHE["x"] = NHE["x"] - n_pos[0] + N[0]
        NHE["y"] = NHE["y"] - n_pos[1] + N[1]
        NHE["z"] = NHE["z"] - n_pos[2] + N[2]
        NHE["resid"] += last_resid

        df = pd.concat([df, NHE], ignore_index=True)
        df["anum"] = range(1, len(df) + 1)
        ####### Rotate NME group #######
        CA = (
            df[(df["atom_name"] == "CA") & (df["resid"] == last_resid)]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        HN1 = (
            df[(df["aa_name"] == "NHE") & (df["atom_name"] == "HN1")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        HN2 = (
            df[(df["aa_name"] == "NHE") & (df["atom_name"] == "HN2")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )
        O = (
            df[(df["resid"] == last_resid) & (df["atom_name"] == "O")]
            .iloc[0][["x", "y", "z"]]
            .to_numpy()
        )

        C = np.asarray(C, dtype=np.float64)
        CA = np.asarray(CA, dtype=np.float64)
        N = np.asarray(N, dtype=np.float64)
        HN1 = np.asarray(HN1, dtype=np.float64)
        HN2 = np.asarray(HN2, dtype=np.float64)
        O = np.asarray(O, dtype=np.float64)

        v1 = C - N
        v2 = HN1 - N
        v3 = rotate_vector_plane_3D(v1, v2, np.pi * 2 / 3) * np.linalg.norm(
            v2
        )  # rotate v2 to v1 by 120 degree
        v4 = rotate_vector_plane_3D(v1, v2, np.pi * 4 / 3) * np.linalg.norm(
            HN2 - N
        )  # also get position of H

        HN1 = N + v3  # correct position of CH3, N
        HN2 = N + v4

        df.loc[
            (df["atom_name"] == "HN1") & (df["aa_name"] == "NHE"), ["x", "y", "z"]
        ] = HN1
        df.loc[
            (df["atom_name"] == "HN2") & (df["aa_name"] == "NHE"), ["x", "y", "z"]
        ] = HN2

        # ################## calculate H coordinates #########
        # p1 = H
        # p2 = N
        # p3 = CH3

        # # Parameters
        # dihedral_angle_deg = 97.150; bond_angle_deg = 110.00; bond_length = 1.078

        # h_coordinate = place_fourth_atom(p1, p2, p3, bond_length, bond_angle_deg,
        # dihedral_angle_deg) ## important

        # angles = [120, 240]  # Angles to rotate by
        # rotated_h_coords = rotate_around_axis(p2, p3, h_coordinate, angles)

        # df.loc[(df['atom_name'] == 'HH31') & (df['aa_name'] == 'NME'), ['x', 'y',
        # 'z']] = h_coordinate
        # df.loc[(df['atom_name'] == 'HH32') & (df['aa_name'] == 'NME'), ['x', 'y',
        # 'z']] = rotated_h_coords[0]
        # df.loc[(df['atom_name'] == 'HH33') & (df['aa_name'] == 'NME'), ['x', 'y',
        # 'z']] = rotated_h_coords[1]

        # #################### rotate C(amino acid) - N(NME) axis to 180 degree
        # ################
        theta = dihedral_angle(O, C, N, HN1)
        rotation_angle = 0 * np.pi - theta

        rotation_axis = N - C
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize axis
        rotation_center = N

        # Copy NME group and rotate its atoms
        NHE_atoms = df[df["aa_name"] == "NHE"].copy()

        for i, row in NHE_atoms.iterrows():
            if row["atom_name"] != "N":
                v = (
                    np.array(row[["x", "y", "z"]].values, dtype=np.float64)
                    - rotation_center
                )  # Translate to origin
                v_rot = rotate_vector(v, rotation_axis, rotation_angle)  # Rotate
                NHE_atoms.loc[i, ["x", "y", "z"]] = (
                    v_rot + rotation_center
                ).tolist()  # Translate back

        # Replace old NME with rotated version
        df = df[df["aa_name"] != "NHE"]
        df = pd.concat([df, NHE_atoms], ignore_index=True)

        # Keep terminal O fixed; the C-terminal cap N replaces the removed OXT.

    # ----------------------------------------------------------------------------------

    ############## correct ACE position ##############
    CA = (
        df[(df["atom_name"] == "CA") & (df["resid"] == 2)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )
    N = (
        df[(df["atom_name"] == "N") & (df["resid"] == 2)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )
    H = (
        df[(df["atom_name"] == "H") & (df["resid"] == 2)]
        .iloc[0][["x", "y", "z"]]
        .to_numpy()
    )
    CA = np.asarray(CA, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64)
    H = np.asarray(H, dtype=np.float64)

    v1 = H - N
    v2 = CA - N
    v3 = (
        rotate_vector_plane_3D(v1, v2, np.pi * 2 / 3) * 1.290
    )  # rotate v2 to v1 by 120 degree
    C = N + v3  # correct position of CH3, N

    df.loc[(df["atom_name"] == "C") & (df["aa_name"] == "ACE"), ["x", "y", "z"]] = C

    v1 = N - C
    # Use the CA-N direction as the ACE reference; H-N swaps the ACE O/CH3 sides.
    v2 = CA - N
    v3 = rotate_vector_plane_3D(v1, v2, np.pi * 4 / 3) * 1.238
    v4 = rotate_vector_plane_3D(v1, v2, np.pi * 2 / 3) * 1.514
    O = C + v3
    CH3 = C + v4
    df.loc[(df["atom_name"] == "O") & (df["aa_name"] == "ACE"), ["x", "y", "z"]] = O
    df.loc[(df["atom_name"] == "CH3") & (df["aa_name"] == "ACE"), ["x", "y", "z"]] = CH3

    ####### place H3 #########
    p1 = O
    p2 = C
    p3 = CH3

    # Parameters
    dihedral_angle_deg = 60
    bond_angle_deg = 110.00
    bond_length = 1.080

    h_coordinate = place_fourth_atom(
        p1, p2, p3, bond_length, bond_angle_deg, dihedral_angle_deg
    )

    angles = [120, 240]  # Angles to rotate by
    rotated_h_coords = rotate_around_axis(p2, p3, h_coordinate, angles)

    df.loc[(df["atom_name"] == "H1") & (df["aa_name"] == "ACE"), ["x", "y", "z"]] = (
        h_coordinate
    )
    df.loc[(df["atom_name"] == "H2") & (df["aa_name"] == "ACE"), ["x", "y", "z"]] = (
        rotated_h_coords[0]
    )
    df.loc[(df["atom_name"] == "H3") & (df["aa_name"] == "ACE"), ["x", "y", "z"]] = (
        rotated_h_coords[1]
    )

    df["anum"] = range(1, len(df) + 1)
    return df.reset_index(drop=True)


# %% Cell cd6c7834
def rotate_vector(v, axis, angle):
    """Function:
        Rotate one 3D vector around a defined 3D axis by a defined angle.

    Parameters
    ----------
    v : numpy.ndarray with shape (3,)
        Vector to rotate.
    axis : numpy.ndarray with shape (3,)
        Rotation axis, for example x=[1,0,0], y=[0,1,0], or z=[0,0,1].
    angle : float
        Rotation angle in radians.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Rotated vector.
    """
    axis = axis / np.linalg.norm(axis)  # Normalize the axis
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    cross_prod = np.cross(axis, v)
    dot_prod = np.dot(axis, v)
    return v * cos_theta + cross_prod * sin_theta + axis * dot_prod * (1 - cos_theta)


def rotate_around_axis(p1, p2, p3, angles):
    """Function:
        Generate several positions of one point by rotating it around the axis from p1
        to p2.

    Parameters
    ----------
    p1 : numpy.ndarray with shape (3,)
        Coordinate of reference point 1.
    p2 : numpy.ndarray with shape (3,)
        Coordinate of reference point 2.
    p3 : numpy.ndarray with shape (3,)
        Coordinate of reference point 3.
    angles : list or numpy.ndarray of floats in degrees
        One output position is made for each angle.

    Returns
    -------
    numpy.ndarray
        Array with shape (len(angles), 3).
    """
    # Calculate the p1-p2 axis
    axis = p2 - p1

    # Translate coordinates to make p1 the origin
    p3_translated = p3 - p1

    # Perform rotations and translate back
    rotated_p3_coords = []
    for angle in angles:
        angle_rad = np.radians(angle)  # Convert degrees to radians
        rotated_p3 = rotate_vector(p3_translated, axis, angle_rad)
        rotated_p3_coords.append(rotated_p3 + p1)  # Translate back

    return np.array(rotated_p3_coords)


def rotate_around_axis_one_angle(p1, p2, points, angle):
    """Function:
        Rotate a group of 3D points around the axis from p1 to p2 by one angle.

    Parameters
    ----------
    p1 : numpy.ndarray with shape (3,)
        Coordinate of reference point 1.
    p2 : numpy.ndarray with shape (3,)
        Coordinate of reference point 2.
    points : numpy.ndarray with shape (n, 3)
        Each row is one 3D point.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    numpy.ndarray with shape (n, 3)
        Rotated points.
    """
    axis = p2 - p1  # Define rotation axis
    axis = np.asarray(axis, dtype=np.float64)
    points_translated = points - p1  # Translate points to make p1 the origin
    rotated_points = []

    for p in points_translated:

        p = np.asarray(p, dtype=np.float64)
        rp = rotate_vector(p, axis, angle)
        rotated_points.append(rp + p1)

    return np.array(rotated_points)


def rotate_vector_plane_3D(v, u, angle):
    """Function:
        Rotate and normalize vector v in the plane defined by vectors u and v.

    Parameters
    ----------
    v : numpy.ndarray with shape (3,)
        Vector to rotate.
    u : numpy.ndarray with shape (3,)
        Second vector used with v to define the rotation plane.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Normalized rotated vector.
    """
    n = np.cross(u, v)
    n = np.asarray(n, dtype=np.float64)

    n /= np.linalg.norm(n)  # Normalize
    v_rot = (
        v * np.cos(angle)
        + np.cross(n, v) * np.sin(angle)
        + n * np.dot(n, v) * (1 - np.cos(angle))
    )
    v_rot /= np.linalg.norm(v_rot)

    return v_rot


def dihedral_angle(p1, p2, p3, p4):
    """Function:
        Calculate the signed dihedral angle formed by four 3D points.

    Parameters
    ----------
    p1 : numpy.ndarray with shape (3,)
        Coordinate of reference point 1.
    p2 : numpy.ndarray with shape (3,)
        Coordinate of reference point 2.
    p3 : numpy.ndarray with shape (3,)
        Coordinate of reference point 3.
    p4 : numpy.ndarray with shape (3,)
        Coordinate of reference point 4.

    Returns
    -------
    float
        Signed dihedral angle in radians.
    """

    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    cos_phi = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
    sin_phi = np.dot(np.cross(n1, n2), b2) / (
        np.linalg.norm(b2) * np.linalg.norm(n1) * np.linalg.norm(n2)
    )

    phi = np.arctan2(sin_phi, cos_phi)

    return phi


def place_fourth_atom(p1, p2, p3, bond_length, bond_angle_deg, dihedral_deg):
    """Function:
        Place a new atom from three reference atoms, one bond length, one bond angle,
        and one dihedral angle.

    Parameters
    ----------
    p1 : numpy.ndarray with shape (3,)
        Coordinate of reference point 1.
    p2 : numpy.ndarray with shape (3,)
        Coordinate of reference point 2.
    p3 : numpy.ndarray with shape (3,)
        Coordinate of reference point 3.
    bond_length : float
        Bond length in Angstrom.
    bond_angle_deg : float
        Bond angle in degrees.
    dihedral_deg : float
        Dihedral angle in degrees.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Coordinate of the fourth atom.
    """
    # Convert to numpy arrays
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    # Convert angles to radians
    bond_angle = np.deg2rad(bond_angle_deg)
    dihedral = np.deg2rad(dihedral_deg)

    # Build local coordinate system at p3 (N)
    v12 = p1 - p2
    v23 = p2 - p3
    # Normalize v23
    e1 = v23 / np.linalg.norm(v23)

    # Normal to the plane
    n = np.cross(v12, v23)
    e2 = n / np.linalg.norm(n)

    # Orthogonal to both e1 and e2
    e3 = np.cross(e1, e2)

    # Coordinates in the local frame
    x = bond_length * np.cos(bond_angle)
    y = -bond_length * np.sin(bond_angle) * np.cos(dihedral + np.pi / 2)
    z = bond_length * np.sin(bond_angle) * np.sin(dihedral + np.pi / 2)

    # Position of H in global coordinates
    p4 = p3 + x * e1 + y * e2 + z * e3
    return p4


# %% Cell local-add-hydrogen
def _read_simple_pdb(pdbfile):
    """Function:
        Read ATOM records from a PDB file into the atom-table format used in this
        notebook.

    Parameters
    ----------
    pdbfile : str or pathlib.Path
        Input PDB file.

    Returns
    -------
    pandas.DataFrame
        One row per ATOM record with atom, anum, atom_name, aa_name, resid, x, y, z.
    """
    records = []
    with open(pdbfile, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            records.append(
                {
                    "atom": line[0:6].strip() or "ATOM",
                    "anum": int(line[6:11]),
                    "atom_name": line[11:17].strip(),
                    "aa_name": line[17:21].strip(),
                    "resid": int(line[22:26]),
                    "x": float(line[26:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            )
    return pd.DataFrame(
        records,
        columns=["atom", "anum", "atom_name", "aa_name", "resid", "x", "y", "z"],
    )


def _place_opposite_bisector(center, neighbors, bond_length):
    """Function:
        Place a new atom opposite to the average direction of its bonded neighboring
        atoms.

    Parameters
    ----------
    center : numpy.ndarray with shape (3,)
        Xyz center of the operation.
    neighbors : list of numpy.ndarray coordinates
        Atoms bonded to the center atom.
    bond_length : float
        Bond length in Angstrom.

    Returns
    -------
    numpy.ndarray or None
        Array with shape (3,), or None when no direction can be defined.
    """
    direction = np.zeros(3, dtype=np.float64)
    for neighbor in neighbors:
        unit = unit_vector(neighbor - center, allow_zero=True)
        if unit is not None:
            direction += unit
    direction = unit_vector(-direction, allow_zero=True)
    if direction is None:
        return None
    return center + bond_length * direction


# %% Cell local-peptide-builder
# Local PeptideBuilder compatibility layer, vendored from PeptideBuilder 1.1.0.
# Original authors: Matthew Z. Tien, Dariya K. Sydykova, Austin G. Meyer,
# and Claus O. Wilke. The upstream files state that they are provided under
# the MIT License. This cell removes the external PeptideBuilder package
# dependency while preserving the same peptide-building geometry.


class Geo:
    """Geometry base class"""

    residue_name: str

    # Geometry to bring together residue
    peptide_bond: float
    CA_C_N_angle: float
    C_N_CA_angle: float

    # Backbone coordinates
    N_CA_C_angle: float
    CA_N_length: float
    CA_C_length: float
    phi: float
    psi_im1: float
    omega: float

    # Carbonyl atom
    C_O_length: float
    CA_C_O_angle: float
    N_CA_C_O_diangle: float

    def __repr__(self) -> str:
        repr = ""
        for var in self.__dict__:
            repr += "%s = %s\n" % (var, self.__dict__[var])
        return repr


class GlyGeo(Geo):
    """Geometry of Glycine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8914

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5117
        self.N_CA_C_O_diangle = 180.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.residue_name = "G"


class AlaGeo(Geo):
    """Geometry of Alanin"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.068

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.5

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6860

        self.residue_name = "A"


class SerGeo(Geo):
    """Geometry of Serine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.2812

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6618

        self.CB_OG_length = 1.417
        self.CA_CB_OG_angle = 110.773
        self.N_CA_CB_OG_diangle = -63.3

        self.residue_name = "S"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_OG_diangle = rotamers[0]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_OG_diangle = -63.3


class CysGeo(Geo):
    """Geometry of Cystine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8856

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.5037

        self.CB_SG_length = 1.808
        self.CA_CB_SG_angle = 113.8169
        self.N_CA_CB_SG_diangle = -62.2

        self.residue_name = "C"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_SG_diangle = rotamers[0]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_SG_diangle = -62.2


class ValGeo(Geo):
    """Geometry of Valine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 109.7698

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5686
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 123.2347

        self.CB_CG1_length = 1.527
        self.CA_CB_CG1_angle = 110.7
        self.N_CA_CB_CG1_diangle = 177.2

        self.CB_CG2_length = 1.527
        self.CA_CB_CG2_angle = 110.4
        self.N_CA_CB_CG2_diangle = -63.3

        self.residue_name = "V"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG1_diangle = rotamers[0]
            self.N_CA_CB_CG2_diangle = rotamers[1]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG1_diangle = 177.2
            self.N_CA_CB_CG2_dianlge = -63.3


class IleGeo(Geo):
    """Geometry of Isoleucine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 109.7202

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5403
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 123.2347

        self.CB_CG1_length = 1.527
        self.CA_CB_CG1_angle = 110.7
        self.N_CA_CB_CG1_diangle = 59.7

        self.CB_CG2_length = 1.527
        self.CA_CB_CG2_angle = 110.4
        self.N_CA_CB_CG2_diangle = -61.6

        self.CG1_CD1_length = 1.52
        self.CB_CG1_CD1_angle = 113.97
        self.CA_CB_CG1_CD1_diangle = 169.8

        self.residue_name = "I"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG1_diangle = rotamers[0]
            self.N_CA_CB_CG2_diangle = rotamers[1]
            self.CA_CB_CG1_CD1_diangle = rotamers[2]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG1_diangle = -61.6
            self.N_CA_CB_CG2_diangle = 59.7
            self.CA_CB_CG1_CD1_diangle = 169.8


# The following function is commented out, because it is not
# recommended to randomize rotamers for isoleucine. The underlying
# reason for this recommendation is that isoleucine's beta-carbon
# is a chiral center.
##    def generateRandomRotamers(self):
##        rotamer_bins = [-60, 60, 180]
##        tempList = []
##        for i in range(0, 3):
##            tempList.append(random.choice(rotamer_bins))
##        self.inputRotamers(tempList)


class LeuGeo(Geo):
    """Geometry of Leucine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8652

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4647
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.4948

        self.CB_CG_length = 1.53
        self.CA_CB_CG_angle = 116.10
        self.N_CA_CB_CG_diangle = -60.1

        self.CG_CD1_length = 1.524
        self.CB_CG_CD1_angle = 110.27
        self.CA_CB_CG_CD1_diangle = 174.9

        self.CG_CD2_length = 1.525
        self.CB_CG_CD2_angle = 110.58
        self.CA_CB_CG_CD2_diangle = 66.7

        self.residue_name = "L"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD1_diangle = rotamers[1]
            self.CA_CB_CG_CD2_diangle = rotamers[2]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -60.1
            self.CA_CB_CG_CD1_diangle = 174.9
            self.CA_CB_CG_CD2_diangle = 66.7

    def generateRandomRotamers(self):
        rotamer_bins = [-60, 60, 180]
        tempList = []
        for i in range(0, 3):
            tempList.append(random.choice(rotamer_bins))
        self.inputRotamers(tempList)


class ThrGeo(Geo):
    """Geometry of Threonine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.7014

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5359
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 123.0953

        self.CB_OG1_length = 1.43
        self.CA_CB_OG1_angle = 109.18
        self.N_CA_CB_OG1_diangle = 60.0

        self.CB_CG2_length = 1.53
        self.CA_CB_CG2_angle = 111.13
        self.N_CA_CB_CG2_diangle = -60.3

        self.residue_name = "T"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_OG1_diangle = rotamers[0]
            self.N_CA_CB_OG2_diangle = rotamers[1]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_OG1_diangle = -60.3
            self.N_CA_CB_OG2_diangle = 60.0


class ArgGeo(Geo):
    """Geometry of Arginine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.98

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.54
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.76

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 113.83
        self.N_CA_CB_CG_diangle = -65.2

        self.CG_CD_length = 1.52
        self.CB_CG_CD_angle = 111.79
        self.CA_CB_CG_CD_diangle = -179.2

        self.CD_NE_length = 1.46
        self.CG_CD_NE_angle = 111.68
        self.CB_CG_CD_NE_diangle = -179.3

        self.NE_CZ_length = 1.33
        self.CD_NE_CZ_angle = 124.79
        self.CG_CD_NE_CZ_diangle = -178.7

        self.CZ_NH1_length = 1.33
        self.NE_CZ_NH1_angle = 120.64
        self.CD_NE_CZ_NH1_diangle = 0.0

        self.CZ_NH2_length = 1.33
        self.NE_CZ_NH2_angle = 119.63
        self.CD_NE_CZ_NH2_diangle = 180.0

        self.residue_name = "R"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD_diangle = rotamers[1]
            self.CB_CG_CD_NE_diangle = rotamers[2]
            self.CG_CD_NE_CZ_diangle = rotamers[3]
            self.CD_NE_CZ_NH1_diangle = rotamers[4]
            self.CD_NE_CZ_NH2_diangle = rotamers[5]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -65.2
            self.CA_CB_CG_CD_diangle = -179.2
            self.CB_CG_CD_NE_diangle = -179.3
            self.CG_CD_NE_CZ_diangle = -178.7
            self.CD_NE_CZ_NH1_diangle = 0.0
            self.CD_NE_CZ_NH2_diangle = 180.0

    def generateRandomRotamers(self):
        rotamer_bins = [-60, 60, 180]
        tempList = []
        for i in range(0, 6):
            tempList.append(random.choice(rotamer_bins))
        self.inputRotamers(tempList)


class LysGeo(Geo):
    """Geometry of Lysine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.08

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.54
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.76

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 113.83
        self.N_CA_CB_CG_diangle = -64.5

        self.CG_CD_length = 1.52
        self.CB_CG_CD_angle = 111.79
        self.CA_CB_CG_CD_diangle = -178.1

        self.CD_CE_length = 1.46
        self.CG_CD_CE_angle = 111.68
        self.CB_CG_CD_CE_diangle = -179.6

        self.CE_NZ_length = 1.33
        self.CD_CE_NZ_angle = 124.79
        self.CG_CD_CE_NZ_diangle = 179.6

        self.residue_name = "K"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD_diangle = rotamers[1]
            self.CB_CG_CD_CE_diangle = rotamers[2]
            self.CG_CD_CE_NZ_diangle = rotamers[3]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -64.5
            self.CA_CB_CG_CD_diangle = -178.1
            self.CB_CG_CD_CE_diangle = -179.6
            self.CG_CD_CE_NZ_diangle = 179.6

    def generateRandomRotamers(self):
        rotamer_bins = [-60, 60, 180]
        tempList = []
        for i in range(0, 4):
            tempList.append(random.choice(rotamer_bins))
        self.inputRotamers(tempList)


class AspGeo(Geo):
    """Geometry of Aspartic Acid"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.03

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.51
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.82

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 113.06
        self.N_CA_CB_CG_diangle = -66.4

        self.CG_OD1_length = 1.25
        self.CB_CG_OD1_angle = 119.22
        self.CA_CB_CG_OD1_diangle = -46.7

        self.CG_OD2_length = 1.25
        self.CB_CG_OD2_angle = 118.218
        self.CA_CB_CG_OD2_diangle = 180 + self.CA_CB_CG_OD1_diangle

        self.residue_name = "D"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_OD1_diangle = rotamers[1]
            if self.CA_CB_CG_OD1_diangle > 0:
                self.CA_CB_CG_OD2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_OD2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -66.4
            self.CA_CB_CG_OD1_diangle = -46.7
            self.CA_CB_CG_OD2_diangle = 180 + self.CA_CB_CG_OD1_diangle


class AsnGeo(Geo):
    """Geometry of Asparagine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.5

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4826
        self.N_CA_C_O_diangle = -60.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 123.2254

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 112.62
        self.N_CA_CB_CG_diangle = -65.5

        self.CG_OD1_length = 1.23
        self.CB_CG_OD1_angle = 120.85
        self.CA_CB_CG_OD1_diangle = -58.3

        self.CG_ND2_length = 1.33
        self.CB_CG_ND2_angle = 116.48
        self.CA_CB_CG_ND2_diangle = 180.0 + self.CA_CB_CG_OD1_diangle

        self.residue_name = "N"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_OD1_diangle = rotamers[1]
            if self.CA_CB_CG_OD1_diangle > 0:
                self.CA_CB_CG_ND2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_ND2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -65.5
            self.CA_CB_CG_OD1_diangle = -58.3
            self.CA_CB_CG_ND2_diangle = 180.0 + self.CA_CB_CG_OD1_diangle


class GluGeo(Geo):
    """Geometry of Glutamic Acid"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.1703

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.511
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.8702

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 113.82
        self.N_CA_CB_CG_diangle = -63.8

        self.CG_CD_length = 1.52
        self.CB_CG_CD_angle = 113.31
        self.CA_CB_CG_CD_diangle = -179.8

        self.CD_OE1_length = 1.25
        self.CG_CD_OE1_angle = 119.02
        self.CB_CG_CD_OE1_diangle = -6.2

        self.CD_OE2_length = 1.25
        self.CG_CD_OE2_angle = 118.08
        self.CB_CG_CD_OE2_diangle = 180.0 + self.CB_CG_CD_OE1_diangle

        self.residue_name = "E"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD_diangle = rotamers[1]
            self.CB_CG_CD_OE1_diangle = rotamers[2]
            if self.CB_CG_CD_OE1_diangle > 0:
                self.CB_CG_CD_OE2_diangle = rotamers[2] - 180.0
            else:
                self.CB_CG_CD_OE2_diangle = rotamers[2] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -63.8
            self.CA_CB_CG_CD_diangle = -179.8
            self.CB_CG_CD_OE1_diangle = -6.2
            self.CB_CG_CD_OE2_diangle = 180.0 + self.CB_CG_CD_OE1_diangle

    def generateRandomRotamers(self):
        rotamer_bins = [-60, 60, 180]
        tempList = []
        for i in range(0, 3):
            tempList.append(random.choice(rotamer_bins))
        self.inputRotamers(tempList)


class GlnGeo(Geo):
    """Geometry of Glutamine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.0849

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5029
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.8134

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 113.75
        self.N_CA_CB_CG_diangle = -60.2

        self.CG_CD_length = 1.52
        self.CB_CG_CD_angle = 112.78
        self.CA_CB_CG_CD_diangle = -69.6

        self.CD_OE1_length = 1.24
        self.CG_CD_OE1_angle = 120.86
        self.CB_CG_CD_OE1_diangle = -50.5

        self.CD_NE2_length = 1.33
        self.CG_CD_NE2_angle = 116.50
        self.CB_CG_CD_NE2_diangle = 180 + self.CB_CG_CD_OE1_diangle

        self.residue_name = "Q"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD_diangle = rotamers[1]
            self.CB_CG_CD_OE1_diangle = rotamers[2]
            if self.CB_CG_CD_OE1_diangle > 0:
                self.CB_CG_CD_NE2_diangle = rotamers[2] - 180.0
            else:
                self.CB_CG_CD_NE2_diangle = rotamers[2] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -60.2
            self.CA_CB_CG_CD_diangle = -69.6
            self.CB_CG_CD_OE1_diangle = -50.5
            self.CB_CG_CD_NE2_diangle = 180 + self.CB_CG_CD_OE1_diangle

    def generateRandomRotamers(self):
        rotamer_bins = [-60, 60, 180]
        tempList = []
        for i in range(0, 3):
            tempList.append(random.choice(rotamer_bins))
        self.inputRotamers(tempList)


class MetGeo(Geo):
    """Geometry of Methionine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.9416

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4816
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6733

        self.CB_CG_length = 1.52
        self.CA_CB_CG_angle = 113.68
        self.N_CA_CB_CG_diangle = -64.4

        self.CG_SD_length = 1.81
        self.CB_CG_SD_angle = 112.69
        self.CA_CB_CG_SD_diangle = -179.6

        self.SD_CE_length = 1.79
        self.CG_SD_CE_angle = 100.61
        self.CB_CG_SD_CE_diangle = 70.1

        self.residue_name = "M"

    def inputRotamers(self, rotamer: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamer[0]
            self.CA_CB_CG_SD_diangle = rotamer[1]
            self.CB_CG_SD_CE_diangle = rotamer[2]
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -64.4
            self.CA_CB_CG_SD_diangle = -179.6
            self.CB_CG_SD_CE_diangle = 70.1

    def generateRandomRotamers(self):
        rotamer_bins = [-60, 60, 180]
        tempList = []
        for i in range(0, 3):
            tempList.append(random.choice(rotamer_bins))
        self.inputRotamers(tempList)


class HisGeo(Geo):
    """Geometry of Histidine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.0859

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.4732
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6711

        self.CB_CG_length = 1.49
        self.CA_CB_CG_angle = 113.74
        self.N_CA_CB_CG_diangle = -63.2

        self.CG_ND1_length = 1.38
        self.CB_CG_ND1_angle = 122.85
        self.CA_CB_CG_ND1_diangle = -75.7

        self.CG_CD2_length = 1.35
        self.CB_CG_CD2_angle = 130.61
        self.CA_CB_CG_CD2_diangle = 180.0 + self.CA_CB_CG_ND1_diangle

        self.ND1_CE1_length = 1.32
        self.CG_ND1_CE1_angle = 108.5
        self.CB_CG_ND1_CE1_diangle = 180.0

        self.CD2_NE2_length = 1.35
        self.CG_CD2_NE2_angle = 108.5
        self.CB_CG_CD2_NE2_diangle = 180.0

        self.residue_name = "H"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_ND1_diangle = rotamers[1]
            if self.CA_CB_CG_ND1_diangle > 0:
                self.CA_CB_CG_CD2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_CD2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -63.2
            self.CA_CB_CG_ND1_diangle = -75.7
            self.CA_CB_CG_CD2_diangle = 180.0 + self.CA_CB_CG_ND1_diangle


class ProGeo(Geo):
    """Geometry of Proline"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 112.7499

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.2945
        self.N_CA_C_O_diangle = -45.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 115.2975

        self.CB_CG_length = 1.49
        self.CA_CB_CG_angle = 104.21
        self.N_CA_CB_CG_diangle = 29.6

        self.CG_CD_length = 1.50
        self.CB_CG_CD_angle = 105.03
        self.CA_CB_CG_CD_diangle = -34.8

        self.residue_name = "P"


class PheGeo(Geo):
    """Geometry of Phenylalanine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.7528

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5316
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6054

        self.CB_CG_length = 1.50
        self.CA_CB_CG_angle = 113.85
        self.N_CA_CB_CG_diangle = -64.7

        self.CG_CD1_length = 1.39
        self.CB_CG_CD1_angle = 120.0
        self.CA_CB_CG_CD1_diangle = 93.3

        self.CG_CD2_length = 1.39
        self.CB_CG_CD2_angle = 120.0
        self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle - 180.0

        self.CD1_CE1_length = 1.39
        self.CG_CD1_CE1_angle = 120.0
        self.CB_CG_CD1_CE1_diangle = 180.0

        self.CD2_CE2_length = 1.39
        self.CG_CD2_CE2_angle = 120.0
        self.CB_CG_CD2_CE2_diangle = 180.0

        self.CE1_CZ_length = 1.39
        self.CD1_CE1_CZ_angle = 120.0
        self.CG_CD1_CE1_CZ_diangle = 0.0

        self.residue_name = "F"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD1_diangle = rotamers[1]
            if self.CA_CB_CG_CD1_diangle > 0:
                self.CA_CB_CG_CD2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_CD2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -64.7
            self.CA_CB_CG_CD1_diangle = 93.3
            self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle - 180.0


class TyrGeo(Geo):
    """Geometry of Tyrosine"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.9288

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5434
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6023

        self.CB_CG_length = 1.51
        self.CA_CB_CG_angle = 113.8
        self.N_CA_CB_CG_diangle = -64.3

        self.CG_CD1_length = 1.39
        self.CB_CG_CD1_angle = 120.98
        self.CA_CB_CG_CD1_diangle = 93.1

        self.CG_CD2_length = 1.39
        self.CB_CG_CD2_angle = 120.82
        self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle + 180.0

        self.CD1_CE1_length = 1.39
        self.CG_CD1_CE1_angle = 120.0
        self.CB_CG_CD1_CE1_diangle = 180.0

        self.CD2_CE2_length = 1.39
        self.CG_CD2_CE2_angle = 120.0
        self.CB_CG_CD2_CE2_diangle = 180.0

        self.CE1_CZ_length = 1.39
        self.CD1_CE1_CZ_angle = 120.0
        self.CG_CD1_CE1_CZ_diangle = 0.0

        self.CZ_OH_length = 1.39
        self.CE1_CZ_OH_angle = 119.78
        self.CD1_CE1_CZ_OH_diangle = 180.0

        self.residue_name = "Y"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD1_diangle = rotamers[1]
            if self.CA_CB_CG_CD1_diangle > 0:
                self.CA_CB_CG_CD2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_CD2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -64.3
            self.CA_CB_CG_CD1_diangle = 93.1
            self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle + 180.0


class TrpGeo(Geo):
    """Geometry of Tryptophan"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 110.8914

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5117
        self.N_CA_C_O_diangle = 120.0

        self.phi = -120
        self.psi_im1 = 140
        self.omega = 180.0
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 116.642992978143
        self.C_N_CA_angle = 121.382215820277

        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6112

        self.CB_CG_length = 1.50
        self.CA_CB_CG_angle = 114.10
        self.N_CA_CB_CG_diangle = -66.4

        self.CG_CD1_length = 1.37
        self.CB_CG_CD1_angle = 127.07
        self.CA_CB_CG_CD1_diangle = 96.3

        self.CG_CD2_length = 1.43
        self.CB_CG_CD2_angle = 126.66
        self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle - 180.0

        self.CD1_NE1_length = 1.38
        self.CG_CD1_NE1_angle = 108.5
        self.CB_CG_CD1_NE1_diangle = 180.0

        self.CD2_CE2_length = 1.40
        self.CG_CD2_CE2_angle = 108.5
        self.CB_CG_CD2_CE2_diangle = 180.0

        self.CD2_CE3_length = 1.40
        self.CG_CD2_CE3_angle = 133.83
        self.CB_CG_CD2_CE3_diangle = 0.0

        self.CE2_CZ2_length = 1.40
        self.CD2_CE2_CZ2_angle = 120.0
        self.CG_CD2_CE2_CZ2_diangle = 180.0

        self.CE3_CZ3_length = 1.40
        self.CD2_CE3_CZ3_angle = 120.0
        self.CG_CD2_CE3_CZ3_diangle = 180.0

        self.CZ2_CH2_length = 1.40
        self.CE2_CZ2_CH2_angle = 120.0
        self.CD2_CE2_CZ2_CH2_diangle = 0.0

        self.residue_name = "W"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD1_diangle = rotamers[1]
            if self.CA_CB_CG_CD1_diangle > 0:
                self.CA_CB_CG_CD2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_CD2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -66.4
            self.CA_CB_CG_CD1_diangle = 96.3
            self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle - 180.0


def geometry(AA: str) -> Geo:
    """Create the geometry-parameter object for one amino-acid type."""
    if AA == "G":
        return GlyGeo()
    elif AA == "A":
        return AlaGeo()
    elif AA == "S":
        return SerGeo()
    elif AA == "C":
        return CysGeo()
    elif AA == "V":
        return ValGeo()
    elif AA == "I":
        return IleGeo()
    elif AA == "L":
        return LeuGeo()
    elif AA == "T":
        return ThrGeo()
    elif AA == "R":
        return ArgGeo()
    elif AA == "K":
        return LysGeo()
    elif AA == "D":
        return AspGeo()
    elif AA == "E":
        return GluGeo()
    elif AA == "N":
        return AsnGeo()
    elif AA == "Q":
        return GlnGeo()
    elif AA == "M":
        return MetGeo()
    elif AA == "H":
        return HisGeo()
    elif AA == "P":
        return ProGeo()
    elif AA == "F":
        return PheGeo()
    elif AA == "Y":
        return TyrGeo()
    elif AA == "W":
        return TrpGeo()
    else:
        return GlyGeo()


# ---- PeptideBuilder.py implementation, adapted for a single notebook cell ----



class Vector:
    """Small replacement for the Biopython Vector methods used by PeptideBuilder."""
    def __init__(self, *values):
        if len(values) == 1:
            self.array = np.asarray(values[0], dtype=float)
        else:
            self.array = np.asarray(values, dtype=float)

    def __getitem__(self, index):
        return self.array[index]

    def __add__(self, other):
        return Vector(self.array + _as_array(other))

    def __sub__(self, other):
        return Vector(self.array - _as_array(other))

    def __mul__(self, value):
        return Vector(self.array * value)

    def __rmul__(self, value):
        return self.__mul__(value)

    def __neg__(self):
        return Vector(-self.array)

    def get_array(self):
        return self.array.copy()

    def left_multiply(self, matrix):
        return Vector(np.asarray(matrix, dtype=float) @ self.array)


def _as_array(value):
    """Convert a local Vector, Atom, or other coordinate value into a NumPy array."""
    if isinstance(value, Vector):
        return value.array
    if isinstance(value, Atom):
        return value.coord
    return np.asarray(value, dtype=float)


def rotaxis(theta, vector):
    """Create a 3 by 3 rotation matrix around a defined axis."""
    axis = _as_array(vector)
    norm = np.linalg.norm(axis)
    if norm == 0:
        raise ValueError("Cannot build a rotation matrix around a zero-length axis.")
    axis = axis / norm
    x, y, z = axis
    c = math.cos(theta)
    s = math.sin(theta)
    C = 1.0 - c
    return np.array([
        [x*x*C + c,   x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, y*y*C + c,   y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, z*z*C + c],
    ])


def calc_angle(v1, v2, v3):
    """Calculate the bond angle formed by three atoms, with v2 as the center atom."""
    a = _as_array(v1) - _as_array(v2)
    b = _as_array(v3) - _as_array(v2)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        raise ValueError("Cannot calculate angle with a zero-length vector.")
    return math.acos(np.clip(np.dot(a, b) / denom, -1.0, 1.0))


def calc_dihedral(v1, v2, v3, v4):
    """Calculate the signed torsion angle formed by four atoms."""
    p0 = _as_array(v1)
    p1 = _as_array(v2)
    p2 = _as_array(v3)
    p3 = _as_array(v4)

    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.atan2(y, x)


class Atom:
    def __init__(self, name, coord, bfactor=0.0, occupancy=1.0, altloc=" ", fullname=None, serial_number=0, element=""):
        self.name = name.strip()
        self.fullname = (fullname or name).strip()
        self.coord = np.asarray(coord, dtype=float)
        self.bfactor = bfactor
        self.occupancy = occupancy
        self.altloc = altloc
        self.serial_number = serial_number
        self.element = (element or self.name[0]).strip()

    def get_vector(self):
        return Vector(self.coord)

    def get_coord(self):
        return self.coord.copy()

    def set_coord(self, coord):
        self.coord = np.asarray(coord, dtype=float)

    def get_id(self):
        return self.name


class Residue:
    def __init__(self, residue_id, resname, segid=""):
        self.id = residue_id
        self.resname = resname
        self.segid = segid
        self.child_list = []
        self.child_dict = {}

    def add(self, atom):
        self.child_list.append(atom)
        self.child_dict[atom.name.strip()] = atom

    def __getitem__(self, atom_name):
        return self.child_dict[atom_name.strip()]

    def get_id(self):
        return self.id


class Chain:
    def __init__(self, chain_id):
        self.id = chain_id
        self.child_list = []
        self.child_dict = {}

    def add(self, residue):
        self.child_list.append(residue)
        self.child_dict[residue.get_id()] = residue

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.child_list[key]
        return self.child_dict[key]


class Model:
    def __init__(self, model_id):
        self.id = model_id
        self.child_list = []
        self.child_dict = {}

    def add(self, chain):
        self.child_list.append(chain)
        self.child_dict[chain.id] = chain

    def __getitem__(self, key):
        return self.child_dict[key]


class Structure:
    def __init__(self, structure_id):
        self.id = structure_id
        self.child_list = []
        self.child_dict = {}

    def add(self, model):
        self.child_list.append(model)
        self.child_dict[model.id] = model

    def __getitem__(self, key):
        return self.child_dict[key]


_AA3_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def is_aa(residue):
    """Check whether a local Residue object is one of the 20 standard amino acids."""
    return getattr(residue, "resname", "") in _AA3_NAMES


class PDBIO:
    """Small PDB writer compatible with the subset used in this notebook."""
    def set_structure(self, structure):
        self.structure = structure

    def save(self, filename):
        serial = 1
        with open(filename, "w") as handle:
            for model in self.structure.child_list:
                for chain in model.child_list:
                    for residue in chain.child_list:
                        resid = residue.get_id()[1]
                        for atom in residue.child_list:
                            x, y, z = atom.coord
                            atom_name = atom.name
                            element = atom.element or atom_name[0]
                            handle.write(
                                f"ATOM  {serial:5d} {atom_name:<4} {residue.resname:>3} {chain.id}{resid:4d}    "
                                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n"
                            )
                            serial += 1
            handle.write("END\n")

def calculateCoordinates(
    refA: Residue, refB: Residue, refC: Residue, L: float, ang: float, di: float
) -> np.ndarray:
    """Calculate a new atom coordinate from three reference atoms and internal-coordinate geometry."""
    AV = refA.get_vector()
    BV = refB.get_vector()
    CV = refC.get_vector()

    CA = AV - CV
    CB = BV - CV

    ##CA vector
    AX = CA[0]
    AY = CA[1]
    AZ = CA[2]

    ##CB vector
    BX = CB[0]
    BY = CB[1]
    BZ = CB[2]

    ##Plane Parameters
    A = (AY * BZ) - (AZ * BY)
    B = (AZ * BX) - (AX * BZ)
    G = (AX * BY) - (AY * BX)

    ##Dot Product Constant
    F = math.sqrt(BX * BX + BY * BY + BZ * BZ) * L * math.cos(ang * (math.pi / 180.0))

    ##Constants
    const = math.sqrt(
        math.pow((B * BZ - BY * G), 2)
        * (
            -(F * F) * (A * A + B * B + G * G)
            + (
                B * B * (BX * BX + BZ * BZ)
                + A * A * (BY * BY + BZ * BZ)
                - (2 * A * BX * BZ * G)
                + (BX * BX + BY * BY) * G * G
                - (2 * B * BY) * (A * BX + BZ * G)
            )
            * L
            * L
        )
    )
    denom = (
        (B * B) * (BX * BX + BZ * BZ)
        + (A * A) * (BY * BY + BZ * BZ)
        - (2 * A * BX * BZ * G)
        + (BX * BX + BY * BY) * (G * G)
        - (2 * B * BY) * (A * BX + BZ * G)
    )

    X = (
        (B * B * BX * F) - (A * B * BY * F) + (F * G) * (-A * BZ + BX * G) + const
    ) / denom

    if (B == 0 or BZ == 0) and (BY == 0 or G == 0):
        const1 = math.sqrt(
            G * G * (-A * A * X * X + (B * B + G * G) * (L - X) * (L + X))
        )
        Y = ((-A * B * X) + const1) / (B * B + G * G)
        Z = -(A * G * G * X + B * const1) / (G * (B * B + G * G))
    else:
        Y = (
            (A * A * BY * F) * (B * BZ - BY * G)
            + G * (-F * math.pow(B * BZ - BY * G, 2) + BX * const)
            - A * (B * B * BX * BZ * F - B * BX * BY * F * G + BZ * const)
        ) / ((B * BZ - BY * G) * denom)
        Z = (
            (A * A * BZ * F) * (B * BZ - BY * G)
            + (B * F) * math.pow(B * BZ - BY * G, 2)
            + (A * BX * F * G) * (-B * BZ + BY * G)
            - B * BX * const
            + A * BY * const
        ) / ((B * BZ - BY * G) * denom)

    # Get the new Vector from the origin
    D = Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp = calc_dihedral(AV, BV, CV, D) * (180.0 / math.pi)

    di = di - temp
    rot = rotaxis(math.pi * (di / 180.0), CV - BV)
    D = (D - BV).left_multiply(rot) + BV

    return D.get_array()


def makeGly(segID: int, N, CA, C, O, geo: Geo) -> Residue:
    """Build one unhydrogenated glycine Residue from backbone atoms and its geometry settings."""
    res = Residue((" ", segID, " "), "GLY", "    ")

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    return res


def makeAla(segID: int, N, CA, C, O, geo: AlaGeo) -> Residue:
    """Build one unhydrogenated alanine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")

    res = Residue((" ", segID, " "), "ALA", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    return res


def makeSer(segID: int, N, CA, C, O, geo: SerGeo) -> Residue:
    """Build one unhydrogenated serine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_OG_length = geo.CB_OG_length
    CA_CB_OG_angle = geo.CA_CB_OG_angle
    N_CA_CB_OG_diangle = geo.N_CA_CB_OG_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    oxygen_g = calculateCoordinates(
        N, CA, CB, CB_OG_length, CA_CB_OG_angle, N_CA_CB_OG_diangle
    )
    OG = Atom("OG", oxygen_g, 0.0, 1.0, " ", " OG", 0, "O")

    ##Create Reside Data Structure
    res = Residue((" ", segID, " "), "SER", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(OG)
    return res


def makeCys(segID: int, N, CA, C, O, geo: CysGeo) -> Residue:
    """Build one unhydrogenated cysteine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_SG_length = geo.CB_SG_length
    CA_CB_SG_angle = geo.CA_CB_SG_angle
    N_CA_CB_SG_diangle = geo.N_CA_CB_SG_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    sulfur_g = calculateCoordinates(
        N, CA, CB, CB_SG_length, CA_CB_SG_angle, N_CA_CB_SG_diangle
    )
    SG = Atom("SG", sulfur_g, 0.0, 1.0, " ", " SG", 0, "S")

    res = Residue((" ", segID, " "), "CYS", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(SG)
    return res


def makeVal(segID: int, N, CA, C, O, geo: ValGeo) -> Residue:
    """Build one unhydrogenated valine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG1_length = geo.CB_CG1_length
    CA_CB_CG1_angle = geo.CA_CB_CG1_angle
    N_CA_CB_CG1_diangle = geo.N_CA_CB_CG1_diangle

    CB_CG2_length = geo.CB_CG2_length
    CA_CB_CG2_angle = geo.CA_CB_CG2_angle
    N_CA_CB_CG2_diangle = geo.N_CA_CB_CG2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g1 = calculateCoordinates(
        N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle
    )
    CG1 = Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    carbon_g2 = calculateCoordinates(
        N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")

    res = Residue((" ", segID, " "), "VAL", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG1)
    res.add(CG2)
    return res


def makeIle(segID: int, N, CA, C, O, geo: IleGeo) -> Residue:
    """Build one unhydrogenated isoleucine Residue from backbone atoms and its geometry settings."""
    ##R-group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG1_length = geo.CB_CG1_length
    CA_CB_CG1_angle = geo.CA_CB_CG1_angle
    N_CA_CB_CG1_diangle = geo.N_CA_CB_CG1_diangle

    CB_CG2_length = geo.CB_CG2_length
    CA_CB_CG2_angle = geo.CA_CB_CG2_angle
    N_CA_CB_CG2_diangle = geo.N_CA_CB_CG2_diangle

    CG1_CD1_length = geo.CG1_CD1_length
    CB_CG1_CD1_angle = geo.CB_CG1_CD1_angle
    CA_CB_CG1_CD1_diangle = geo.CA_CB_CG1_CD1_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g1 = calculateCoordinates(
        N, CA, CB, CB_CG1_length, CA_CB_CG1_angle, N_CA_CB_CG1_diangle
    )
    CG1 = Atom("CG1", carbon_g1, 0.0, 1.0, " ", " CG1", 0, "C")
    carbon_g2 = calculateCoordinates(
        N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")
    carbon_d1 = calculateCoordinates(
        CA, CB, CG1, CG1_CD1_length, CB_CG1_CD1_angle, CA_CB_CG1_CD1_diangle
    )
    CD1 = Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")

    res = Residue((" ", segID, " "), "ILE", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG1)
    res.add(CG2)
    res.add(CD1)
    return res


def makeLeu(segID: int, N, CA, C, O, geo: LeuGeo) -> Residue:
    """Build one unhydrogenated leucine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD1_length = geo.CG_CD1_length
    CB_CG_CD1_angle = geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle = geo.CA_CB_CG_CD1_diangle

    CG_CD2_length = geo.CG_CD2_length
    CB_CG_CD2_angle = geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle = geo.CA_CB_CG_CD2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g1 = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g1, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1 = calculateCoordinates(
        CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle
    )
    CD1 = Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2 = calculateCoordinates(
        CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")

    res = Residue((" ", segID, " "), "LEU", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CD2)
    return res


def makeThr(segID: int, N, CA, C, O, geo: ThrGeo) -> Residue:
    """Build one unhydrogenated threonine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_OG1_length = geo.CB_OG1_length
    CA_CB_OG1_angle = geo.CA_CB_OG1_angle
    N_CA_CB_OG1_diangle = geo.N_CA_CB_OG1_diangle

    CB_CG2_length = geo.CB_CG2_length
    CA_CB_CG2_angle = geo.CA_CB_CG2_angle
    N_CA_CB_CG2_diangle = geo.N_CA_CB_CG2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    oxygen_g1 = calculateCoordinates(
        N, CA, CB, CB_OG1_length, CA_CB_OG1_angle, N_CA_CB_OG1_diangle
    )
    OG1 = Atom("OG1", oxygen_g1, 0.0, 1.0, " ", " OG1", 0, "O")
    carbon_g2 = calculateCoordinates(
        N, CA, CB, CB_CG2_length, CA_CB_CG2_angle, N_CA_CB_CG2_diangle
    )
    CG2 = Atom("CG2", carbon_g2, 0.0, 1.0, " ", " CG2", 0, "C")

    res = Residue((" ", segID, " "), "THR", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(OG1)
    res.add(CG2)
    return res


def makeArg(segID: int, N, CA, C, O, geo: ArgGeo) -> Residue:
    """Build one unhydrogenated arginine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD_length = geo.CG_CD_length
    CB_CG_CD_angle = geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle = geo.CA_CB_CG_CD_diangle

    CD_NE_length = geo.CD_NE_length
    CG_CD_NE_angle = geo.CG_CD_NE_angle
    CB_CG_CD_NE_diangle = geo.CB_CG_CD_NE_diangle

    NE_CZ_length = geo.NE_CZ_length
    CD_NE_CZ_angle = geo.CD_NE_CZ_angle
    CG_CD_NE_CZ_diangle = geo.CG_CD_NE_CZ_diangle

    CZ_NH1_length = geo.CZ_NH1_length
    NE_CZ_NH1_angle = geo.NE_CZ_NH1_angle
    CD_NE_CZ_NH1_diangle = geo.CD_NE_CZ_NH1_diangle

    CZ_NH2_length = geo.CZ_NH2_length
    NE_CZ_NH2_angle = geo.NE_CZ_NH2_angle
    CD_NE_CZ_NH2_diangle = geo.CD_NE_CZ_NH2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d = calculateCoordinates(
        CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    nitrogen_e = calculateCoordinates(
        CB, CG, CD, CD_NE_length, CG_CD_NE_angle, CB_CG_CD_NE_diangle
    )
    NE = Atom("NE", nitrogen_e, 0.0, 1.0, " ", " NE", 0, "N")
    carbon_z = calculateCoordinates(
        CG, CD, NE, NE_CZ_length, CD_NE_CZ_angle, CG_CD_NE_CZ_diangle
    )
    CZ = Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    nitrogen_h1 = calculateCoordinates(
        CD, NE, CZ, CZ_NH1_length, NE_CZ_NH1_angle, CD_NE_CZ_NH1_diangle
    )
    NH1 = Atom("NH1", nitrogen_h1, 0.0, 1.0, " ", " NH1", 0, "N")
    nitrogen_h2 = calculateCoordinates(
        CD, NE, CZ, CZ_NH2_length, NE_CZ_NH2_angle, CD_NE_CZ_NH2_diangle
    )
    NH2 = Atom("NH2", nitrogen_h2, 0.0, 1.0, " ", " NH2", 0, "N")

    res = Residue((" ", segID, " "), "ARG", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(NE)
    res.add(CZ)
    res.add(NH1)
    res.add(NH2)
    return res


def makeLys(segID: int, N, CA, C, O, geo: LysGeo) -> Residue:
    """Build one unhydrogenated lysine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD_length = geo.CG_CD_length
    CB_CG_CD_angle = geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle = geo.CA_CB_CG_CD_diangle

    CD_CE_length = geo.CD_CE_length
    CG_CD_CE_angle = geo.CG_CD_CE_angle
    CB_CG_CD_CE_diangle = geo.CB_CG_CD_CE_diangle

    CE_NZ_length = geo.CE_NZ_length
    CD_CE_NZ_angle = geo.CD_CE_NZ_angle
    CG_CD_CE_NZ_diangle = geo.CG_CD_CE_NZ_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d = calculateCoordinates(
        CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    carbon_e = calculateCoordinates(
        CB, CG, CD, CD_CE_length, CG_CD_CE_angle, CB_CG_CD_CE_diangle
    )
    CE = Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")
    nitrogen_z = calculateCoordinates(
        CG, CD, CE, CE_NZ_length, CD_CE_NZ_angle, CG_CD_CE_NZ_diangle
    )
    NZ = Atom("NZ", nitrogen_z, 0.0, 1.0, " ", " NZ", 0, "N")

    res = Residue((" ", segID, " "), "LYS", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(CE)
    res.add(NZ)
    return res


def makeAsp(segID: int, N, CA, C, O, geo: AspGeo) -> Residue:
    """Build one unhydrogenated aspartic acid Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_OD1_length = geo.CG_OD1_length
    CB_CG_OD1_angle = geo.CB_CG_OD1_angle
    CA_CB_CG_OD1_diangle = geo.CA_CB_CG_OD1_diangle

    CG_OD2_length = geo.CG_OD2_length
    CB_CG_OD2_angle = geo.CB_CG_OD2_angle
    CA_CB_CG_OD2_diangle = geo.CA_CB_CG_OD2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    oxygen_d1 = calculateCoordinates(
        CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle
    )
    OD1 = Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")
    oxygen_d2 = calculateCoordinates(
        CA, CB, CG, CG_OD2_length, CB_CG_OD2_angle, CA_CB_CG_OD2_diangle
    )
    OD2 = Atom("OD2", oxygen_d2, 0.0, 1.0, " ", " OD2", 0, "O")

    res = Residue((" ", segID, " "), "ASP", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(OD1)
    res.add(OD2)
    return res


def makeAsn(segID, N, CA, C, O, geo):
    """Build one unhydrogenated asparagine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_OD1_length = geo.CG_OD1_length
    CB_CG_OD1_angle = geo.CB_CG_OD1_angle
    CA_CB_CG_OD1_diangle = geo.CA_CB_CG_OD1_diangle

    CG_ND2_length = geo.CG_ND2_length
    CB_CG_ND2_angle = geo.CB_CG_ND2_angle
    CA_CB_CG_ND2_diangle = geo.CA_CB_CG_ND2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    oxygen_d1 = calculateCoordinates(
        CA, CB, CG, CG_OD1_length, CB_CG_OD1_angle, CA_CB_CG_OD1_diangle
    )
    OD1 = Atom("OD1", oxygen_d1, 0.0, 1.0, " ", " OD1", 0, "O")
    nitrogen_d2 = calculateCoordinates(
        CA, CB, CG, CG_ND2_length, CB_CG_ND2_angle, CA_CB_CG_ND2_diangle
    )
    ND2 = Atom("ND2", nitrogen_d2, 0.0, 1.0, " ", " ND2", 0, "N")
    res = Residue((" ", segID, " "), "ASN", "    ")

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(OD1)
    res.add(ND2)
    return res


def makeGlu(segID: int, N, CA, C, O, geo: GluGeo) -> Residue:
    """Build one unhydrogenated glutamic acid Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD_length = geo.CG_CD_length
    CB_CG_CD_angle = geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle = geo.CA_CB_CG_CD_diangle

    CD_OE1_length = geo.CD_OE1_length
    CG_CD_OE1_angle = geo.CG_CD_OE1_angle
    CB_CG_CD_OE1_diangle = geo.CB_CG_CD_OE1_diangle

    CD_OE2_length = geo.CD_OE2_length
    CG_CD_OE2_angle = geo.CG_CD_OE2_angle
    CB_CG_CD_OE2_diangle = geo.CB_CG_CD_OE2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d = calculateCoordinates(
        CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    oxygen_e1 = calculateCoordinates(
        CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle
    )
    OE1 = Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    oxygen_e2 = calculateCoordinates(
        CB, CG, CD, CD_OE2_length, CG_CD_OE2_angle, CB_CG_CD_OE2_diangle
    )
    OE2 = Atom("OE2", oxygen_e2, 0.0, 1.0, " ", " OE2", 0, "O")

    res = Residue((" ", segID, " "), "GLU", "    ")

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(OE1)
    res.add(OE2)
    return res


def makeGln(segID: int, N, CA, C, O, geo: GlnGeo) -> Residue:
    """Build one unhydrogenated glutamine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD_length = geo.CG_CD_length
    CB_CG_CD_angle = geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle = geo.CA_CB_CG_CD_diangle

    CD_OE1_length = geo.CD_OE1_length
    CG_CD_OE1_angle = geo.CG_CD_OE1_angle
    CB_CG_CD_OE1_diangle = geo.CB_CG_CD_OE1_diangle

    CD_NE2_length = geo.CD_NE2_length
    CG_CD_NE2_angle = geo.CG_CD_NE2_angle
    CB_CG_CD_NE2_diangle = geo.CB_CG_CD_NE2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d = calculateCoordinates(
        CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")
    oxygen_e1 = calculateCoordinates(
        CB, CG, CD, CD_OE1_length, CG_CD_OE1_angle, CB_CG_CD_OE1_diangle
    )
    OE1 = Atom("OE1", oxygen_e1, 0.0, 1.0, " ", " OE1", 0, "O")
    nitrogen_e2 = calculateCoordinates(
        CB, CG, CD, CD_NE2_length, CG_CD_NE2_angle, CB_CG_CD_NE2_diangle
    )
    NE2 = Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")

    ##Create Residue DS
    res = Residue((" ", segID, " "), "GLN", "    ")

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)
    res.add(OE1)
    res.add(NE2)
    return res


def makeMet(segID: int, N, CA, C, O, geo: MetGeo) -> Residue:
    """Build one unhydrogenated methionine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_SD_length = geo.CG_SD_length
    CB_CG_SD_angle = geo.CB_CG_SD_angle
    CA_CB_CG_SD_diangle = geo.CA_CB_CG_SD_diangle

    SD_CE_length = geo.SD_CE_length
    CG_SD_CE_angle = geo.CG_SD_CE_angle
    CB_CG_SD_CE_diangle = geo.CB_CG_SD_CE_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    sulfur_d = calculateCoordinates(
        CA, CB, CG, CG_SD_length, CB_CG_SD_angle, CA_CB_CG_SD_diangle
    )
    SD = Atom("SD", sulfur_d, 0.0, 1.0, " ", " SD", 0, "S")
    carbon_e = calculateCoordinates(
        CB, CG, SD, SD_CE_length, CG_SD_CE_angle, CB_CG_SD_CE_diangle
    )
    CE = Atom("CE", carbon_e, 0.0, 1.0, " ", " CE", 0, "C")

    res = Residue((" ", segID, " "), "MET", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(SD)
    res.add(CE)
    return res


def makeHis(segID: int, N, CA, C, O, geo: HisGeo) -> Residue:
    """Build one unhydrogenated histidine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_ND1_length = geo.CG_ND1_length
    CB_CG_ND1_angle = geo.CB_CG_ND1_angle
    CA_CB_CG_ND1_diangle = geo.CA_CB_CG_ND1_diangle

    CG_CD2_length = geo.CG_CD2_length
    CB_CG_CD2_angle = geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle = geo.CA_CB_CG_CD2_diangle

    ND1_CE1_length = geo.ND1_CE1_length
    CG_ND1_CE1_angle = geo.CG_ND1_CE1_angle
    CB_CG_ND1_CE1_diangle = geo.CB_CG_ND1_CE1_diangle

    CD2_NE2_length = geo.CD2_NE2_length
    CG_CD2_NE2_angle = geo.CG_CD2_NE2_angle
    CB_CG_CD2_NE2_diangle = geo.CB_CG_CD2_NE2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    nitrogen_d1 = calculateCoordinates(
        CA, CB, CG, CG_ND1_length, CB_CG_ND1_angle, CA_CB_CG_ND1_diangle
    )
    ND1 = Atom("ND1", nitrogen_d1, 0.0, 1.0, " ", " ND1", 0, "N")
    carbon_d2 = calculateCoordinates(
        CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e1 = calculateCoordinates(
        CB, CG, ND1, ND1_CE1_length, CG_ND1_CE1_angle, CB_CG_ND1_CE1_diangle
    )
    CE1 = Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    nitrogen_e2 = calculateCoordinates(
        CB, CG, CD2, CD2_NE2_length, CG_CD2_NE2_angle, CB_CG_CD2_NE2_diangle
    )
    NE2 = Atom("NE2", nitrogen_e2, 0.0, 1.0, " ", " NE2", 0, "N")

    res = Residue((" ", segID, " "), "HIS", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(ND1)
    res.add(CD2)
    res.add(CE1)
    res.add(NE2)
    return res


def makePro(segID: int, N, CA, C, O, geo: ProGeo) -> Residue:
    """Build one unhydrogenated proline Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD_length = geo.CG_CD_length
    CB_CG_CD_angle = geo.CB_CG_CD_angle
    CA_CB_CG_CD_diangle = geo.CA_CB_CG_CD_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d = calculateCoordinates(
        CA, CB, CG, CG_CD_length, CB_CG_CD_angle, CA_CB_CG_CD_diangle
    )
    CD = Atom("CD", carbon_d, 0.0, 1.0, " ", " CD", 0, "C")

    res = Residue((" ", segID, " "), "PRO", "    ")

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD)

    return res


def makePhe(segID: int, N, CA, C, O, geo: PheGeo) -> Residue:
    """Build one unhydrogenated phenylalanine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD1_length = geo.CG_CD1_length
    CB_CG_CD1_angle = geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle = geo.CA_CB_CG_CD1_diangle

    CG_CD2_length = geo.CG_CD2_length
    CB_CG_CD2_angle = geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle = geo.CA_CB_CG_CD2_diangle

    CD1_CE1_length = geo.CD1_CE1_length
    CG_CD1_CE1_angle = geo.CG_CD1_CE1_angle
    CB_CG_CD1_CE1_diangle = geo.CB_CG_CD1_CE1_diangle

    CD2_CE2_length = geo.CD2_CE2_length
    CG_CD2_CE2_angle = geo.CG_CD2_CE2_angle
    CB_CG_CD2_CE2_diangle = geo.CB_CG_CD2_CE2_diangle

    CE1_CZ_length = geo.CE1_CZ_length
    CD1_CE1_CZ_angle = geo.CD1_CE1_CZ_angle
    CG_CD1_CE1_CZ_diangle = geo.CG_CD1_CE1_CZ_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1 = calculateCoordinates(
        CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle
    )
    CD1 = Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2 = calculateCoordinates(
        CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e1 = calculateCoordinates(
        CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle
    )
    CE1 = Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    carbon_e2 = calculateCoordinates(
        CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_z = calculateCoordinates(
        CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle
    )
    CZ = Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")

    res = Residue((" ", segID, " "), "PHE", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CE1)
    res.add(CD2)
    res.add(CE2)
    res.add(CZ)
    return res


def makeTyr(segID: int, N, CA, C, O, geo: TyrGeo) -> Residue:
    """Build one unhydrogenated tyrosine Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD1_length = geo.CG_CD1_length
    CB_CG_CD1_angle = geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle = geo.CA_CB_CG_CD1_diangle

    CG_CD2_length = geo.CG_CD2_length
    CB_CG_CD2_angle = geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle = geo.CA_CB_CG_CD2_diangle

    CD1_CE1_length = geo.CD1_CE1_length
    CG_CD1_CE1_angle = geo.CG_CD1_CE1_angle
    CB_CG_CD1_CE1_diangle = geo.CB_CG_CD1_CE1_diangle

    CD2_CE2_length = geo.CD2_CE2_length
    CG_CD2_CE2_angle = geo.CG_CD2_CE2_angle
    CB_CG_CD2_CE2_diangle = geo.CB_CG_CD2_CE2_diangle

    CE1_CZ_length = geo.CE1_CZ_length
    CD1_CE1_CZ_angle = geo.CD1_CE1_CZ_angle
    CG_CD1_CE1_CZ_diangle = geo.CG_CD1_CE1_CZ_diangle

    CZ_OH_length = geo.CZ_OH_length
    CE1_CZ_OH_angle = geo.CE1_CZ_OH_angle
    CD1_CE1_CZ_OH_diangle = geo.CD1_CE1_CZ_OH_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1 = calculateCoordinates(
        CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle
    )
    CD1 = Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2 = calculateCoordinates(
        CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    carbon_e1 = calculateCoordinates(
        CB, CG, CD1, CD1_CE1_length, CG_CD1_CE1_angle, CB_CG_CD1_CE1_diangle
    )
    CE1 = Atom("CE1", carbon_e1, 0.0, 1.0, " ", " CE1", 0, "C")
    carbon_e2 = calculateCoordinates(
        CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_z = calculateCoordinates(
        CG, CD1, CE1, CE1_CZ_length, CD1_CE1_CZ_angle, CG_CD1_CE1_CZ_diangle
    )
    CZ = Atom("CZ", carbon_z, 0.0, 1.0, " ", " CZ", 0, "C")
    oxygen_h = calculateCoordinates(
        CD1, CE1, CZ, CZ_OH_length, CE1_CZ_OH_angle, CD1_CE1_CZ_OH_diangle
    )
    OH = Atom("OH", oxygen_h, 0.0, 1.0, " ", " OH", 0, "O")

    ##Create Residue Data S
    res = Residue((" ", segID, " "), "TYR", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CE1)
    res.add(CD2)
    res.add(CE2)
    res.add(CZ)
    res.add(OH)
    return res


def makeTrp(segID: int, N, CA, C, O, geo: TrpGeo) -> Residue:
    """Build one unhydrogenated tryptophan Residue from backbone atoms and its geometry settings."""
    ##R-Group
    CA_CB_length = geo.CA_CB_length
    C_CA_CB_angle = geo.C_CA_CB_angle
    N_C_CA_CB_diangle = geo.N_C_CA_CB_diangle

    CB_CG_length = geo.CB_CG_length
    CA_CB_CG_angle = geo.CA_CB_CG_angle
    N_CA_CB_CG_diangle = geo.N_CA_CB_CG_diangle

    CG_CD1_length = geo.CG_CD1_length
    CB_CG_CD1_angle = geo.CB_CG_CD1_angle
    CA_CB_CG_CD1_diangle = geo.CA_CB_CG_CD1_diangle

    CG_CD2_length = geo.CG_CD2_length
    CB_CG_CD2_angle = geo.CB_CG_CD2_angle
    CA_CB_CG_CD2_diangle = geo.CA_CB_CG_CD2_diangle

    CD1_NE1_length = geo.CD1_NE1_length
    CG_CD1_NE1_angle = geo.CG_CD1_NE1_angle
    CB_CG_CD1_NE1_diangle = geo.CB_CG_CD1_NE1_diangle

    CD2_CE2_length = geo.CD2_CE2_length
    CG_CD2_CE2_angle = geo.CG_CD2_CE2_angle
    CB_CG_CD2_CE2_diangle = geo.CB_CG_CD2_CE2_diangle

    CD2_CE3_length = geo.CD2_CE3_length
    CG_CD2_CE3_angle = geo.CG_CD2_CE3_angle
    CB_CG_CD2_CE3_diangle = geo.CB_CG_CD2_CE3_diangle

    CE2_CZ2_length = geo.CE2_CZ2_length
    CD2_CE2_CZ2_angle = geo.CD2_CE2_CZ2_angle
    CG_CD2_CE2_CZ2_diangle = geo.CG_CD2_CE2_CZ2_diangle

    CE3_CZ3_length = geo.CE3_CZ3_length
    CD2_CE3_CZ3_angle = geo.CD2_CE3_CZ3_angle
    CG_CD2_CE3_CZ3_diangle = geo.CG_CD2_CE3_CZ3_diangle

    CZ2_CH2_length = geo.CZ2_CH2_length
    CE2_CZ2_CH2_angle = geo.CE2_CZ2_CH2_angle
    CD2_CE2_CZ2_CH2_diangle = geo.CD2_CE2_CZ2_CH2_diangle

    carbon_b = calculateCoordinates(
        N, C, CA, CA_CB_length, C_CA_CB_angle, N_C_CA_CB_diangle
    )
    CB = Atom("CB", carbon_b, 0.0, 1.0, " ", " CB", 0, "C")
    carbon_g = calculateCoordinates(
        N, CA, CB, CB_CG_length, CA_CB_CG_angle, N_CA_CB_CG_diangle
    )
    CG = Atom("CG", carbon_g, 0.0, 1.0, " ", " CG", 0, "C")
    carbon_d1 = calculateCoordinates(
        CA, CB, CG, CG_CD1_length, CB_CG_CD1_angle, CA_CB_CG_CD1_diangle
    )
    CD1 = Atom("CD1", carbon_d1, 0.0, 1.0, " ", " CD1", 0, "C")
    carbon_d2 = calculateCoordinates(
        CA, CB, CG, CG_CD2_length, CB_CG_CD2_angle, CA_CB_CG_CD2_diangle
    )
    CD2 = Atom("CD2", carbon_d2, 0.0, 1.0, " ", " CD2", 0, "C")
    nitrogen_e1 = calculateCoordinates(
        CB, CG, CD1, CD1_NE1_length, CG_CD1_NE1_angle, CB_CG_CD1_NE1_diangle
    )
    NE1 = Atom("NE1", nitrogen_e1, 0.0, 1.0, " ", " NE1", 0, "N")
    carbon_e2 = calculateCoordinates(
        CB, CG, CD2, CD2_CE2_length, CG_CD2_CE2_angle, CB_CG_CD2_CE2_diangle
    )
    CE2 = Atom("CE2", carbon_e2, 0.0, 1.0, " ", " CE2", 0, "C")
    carbon_e3 = calculateCoordinates(
        CB, CG, CD2, CD2_CE3_length, CG_CD2_CE3_angle, CB_CG_CD2_CE3_diangle
    )
    CE3 = Atom("CE3", carbon_e3, 0.0, 1.0, " ", " CE3", 0, "C")

    carbon_z2 = calculateCoordinates(
        CG, CD2, CE2, CE2_CZ2_length, CD2_CE2_CZ2_angle, CG_CD2_CE2_CZ2_diangle
    )
    CZ2 = Atom("CZ2", carbon_z2, 0.0, 1.0, " ", " CZ2", 0, "C")

    carbon_z3 = calculateCoordinates(
        CG, CD2, CE3, CE3_CZ3_length, CD2_CE3_CZ3_angle, CG_CD2_CE3_CZ3_diangle
    )
    CZ3 = Atom("CZ3", carbon_z3, 0.0, 1.0, " ", " CZ3", 0, "C")

    carbon_h2 = calculateCoordinates(
        CD2, CE2, CZ2, CZ2_CH2_length, CE2_CZ2_CH2_angle, CD2_CE2_CZ2_CH2_diangle
    )
    CH2 = Atom("CH2", carbon_h2, 0.0, 1.0, " ", " CH2", 0, "C")

    ##Create Residue DS
    res = Residue((" ", segID, " "), "TRP", "    ")
    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    res.add(CB)
    res.add(CG)
    res.add(CD1)
    res.add(CD2)

    res.add(NE1)
    res.add(CE2)
    res.add(CE3)

    res.add(CZ2)
    res.add(CZ3)

    res.add(CH2)
    return res


def make_res_of_type(segID: int, N, CA, C, O, geo: Geo) -> Residue:
    """Choose the correct residue-building function from the supplied amino-acid geometry object."""
    if isinstance(geo, GlyGeo):
        res = makeGly(segID, N, CA, C, O, geo)
    elif isinstance(geo, AlaGeo):
        res = makeAla(segID, N, CA, C, O, geo)
    elif isinstance(geo, SerGeo):
        res = makeSer(segID, N, CA, C, O, geo)
    elif isinstance(geo, CysGeo):
        res = makeCys(segID, N, CA, C, O, geo)
    elif isinstance(geo, ValGeo):
        res = makeVal(segID, N, CA, C, O, geo)
    elif isinstance(geo, IleGeo):
        res = makeIle(segID, N, CA, C, O, geo)
    elif isinstance(geo, LeuGeo):
        res = makeLeu(segID, N, CA, C, O, geo)
    elif isinstance(geo, ThrGeo):
        res = makeThr(segID, N, CA, C, O, geo)
    elif isinstance(geo, ArgGeo):
        res = makeArg(segID, N, CA, C, O, geo)
    elif isinstance(geo, LysGeo):
        res = makeLys(segID, N, CA, C, O, geo)
    elif isinstance(geo, AspGeo):
        res = makeAsp(segID, N, CA, C, O, geo)
    elif isinstance(geo, GluGeo):
        res = makeGlu(segID, N, CA, C, O, geo)
    elif isinstance(geo, AsnGeo):
        res = makeAsn(segID, N, CA, C, O, geo)
    elif isinstance(geo, GlnGeo):
        res = makeGln(segID, N, CA, C, O, geo)
    elif isinstance(geo, MetGeo):
        res = makeMet(segID, N, CA, C, O, geo)
    elif isinstance(geo, HisGeo):
        res = makeHis(segID, N, CA, C, O, geo)
    elif isinstance(geo, ProGeo):
        res = makePro(segID, N, CA, C, O, geo)
    elif isinstance(geo, PheGeo):
        res = makePhe(segID, N, CA, C, O, geo)
    elif isinstance(geo, TyrGeo):
        res = makeTyr(segID, N, CA, C, O, geo)
    elif isinstance(geo, TrpGeo):
        res = makeTrp(segID, N, CA, C, O, geo)
    else:
        res = makeGly(segID, N, CA, C, O, geo)

    return res


def initialize_res(residue: Union[Geo, str]) -> Structure:
    """Create a new local peptide Structure containing the first residue."""

    if isinstance(residue, Geo):
        geo = residue
    elif isinstance(residue, str):
        geo = geometry(residue)
    else:
        raise ValueError("Invalid residue argument:", residue)

    segID = 1
    AA = geo.residue_name
    CA_N_length = geo.CA_N_length
    CA_C_length = geo.CA_C_length
    N_CA_C_angle = geo.N_CA_C_angle

    CA_coord = np.array([0.0, 0.0, 0.0])
    C_coord = np.array([CA_C_length, 0, 0])
    N_coord = np.array(
        [
            CA_N_length * math.cos(N_CA_C_angle * (math.pi / 180.0)),
            CA_N_length * math.sin(N_CA_C_angle * (math.pi / 180.0)),
            0,
        ]
    )

    N = Atom("N", N_coord, 0.0, 1.0, " ", " N", 0, "N")
    CA = Atom("CA", CA_coord, 0.0, 1.0, " ", " CA", 0, "C")
    C = Atom("C", C_coord, 0.0, 1.0, " ", " C", 0, "C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length = geo.C_O_length
    CA_C_O_angle = geo.CA_C_O_angle
    N_CA_C_O_diangle = geo.N_CA_C_O_diangle

    carbonyl = calculateCoordinates(
        N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle
    )
    O = Atom("O", carbonyl, 0.0, 1.0, " ", " O", 0, "O")

    res = make_res_of_type(segID, N, CA, C, O, geo)

    cha = Chain("A")
    cha.add(res)

    mod = Model(0)
    mod.add(cha)

    struc = Structure("X")
    struc.add(mod)
    return struc


def getReferenceResidue(structure: Structure) -> Residue:
    """Get the last residue in chain A, which is used as the reference for adding another residue."""

    # If the following line doesn't work we're in trouble.
    # Likely initialize_res() wasn't called.
    resRef = structure[0]["A"].child_list[-1]

    # If the residue is not an amino acid we're in trouble.
    # Likely somebody is trying to append residues to an existing
    # structure that has non-amino-acid molecules in the chain.
    assert is_aa(resRef)

    return resRef


def add_residue_from_geo(structure: Structure, geo: Geo) -> Structure:
    """Add one residue to an existing local Structure using a prepared geometry object."""
    resRef = getReferenceResidue(structure)
    AA = geo.residue_name
    segID = resRef.get_id()[1]
    segID += 1

    ##geometry to bring together residue
    peptide_bond = geo.peptide_bond
    CA_C_N_angle = geo.CA_C_N_angle
    C_N_CA_angle = geo.C_N_CA_angle

    ##Backbone Coordinates
    N_CA_C_angle = geo.N_CA_C_angle
    CA_N_length = geo.CA_N_length
    CA_C_length = geo.CA_C_length
    phi = geo.phi
    psi_im1 = geo.psi_im1
    omega = geo.omega

    N_coord = calculateCoordinates(
        resRef["N"], resRef["CA"], resRef["C"], peptide_bond, CA_C_N_angle, psi_im1
    )
    N = Atom("N", N_coord, 0.0, 1.0, " ", " N", 0, "N")

    CA_coord = calculateCoordinates(
        resRef["CA"], resRef["C"], N, CA_N_length, C_N_CA_angle, omega
    )
    CA = Atom("CA", CA_coord, 0.0, 1.0, " ", " CA", 0, "C")

    C_coord = calculateCoordinates(resRef["C"], N, CA, CA_C_length, N_CA_C_angle, phi)
    C = Atom("C", C_coord, 0.0, 1.0, " ", " C", 0, "C")

    ##Create Carbonyl atom (to be moved later)
    C_O_length = geo.C_O_length
    CA_C_O_angle = geo.CA_C_O_angle
    N_CA_C_O_diangle = geo.N_CA_C_O_diangle

    carbonyl = calculateCoordinates(
        N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle
    )
    O = Atom("O", carbonyl, 0.0, 1.0, " ", " O", 0, "O")

    res = make_res_of_type(segID, N, CA, C, O, geo)

    resRef["O"].set_coord(
        calculateCoordinates(
            res["N"], resRef["CA"], resRef["C"], C_O_length, CA_C_O_angle, 180.0
        )
    )

    ghost = Atom(
        "N",
        calculateCoordinates(
            res["N"], res["CA"], res["C"], peptide_bond, CA_C_N_angle, psi_im1
        ),
        0.0,
        0.0,
        " ",
        "N",
        0,
        "N",
    )
    res["O"].set_coord(
        calculateCoordinates(
            res["N"], res["CA"], res["C"], C_O_length, CA_C_O_angle, 180.0
        )
    )

    structure[0]["A"].add(res)
    return structure


def make_extended_structure(AA_chain: str) -> Structure:
    """Build a peptide Structure from a sequence using the default extended conformation."""
    geo = geometry(AA_chain[0])
    struc = initialize_res(geo)

    for i in range(1, len(AA_chain)):
        AA = AA_chain[i]
        geo = geometry(AA)
        add_residue(struc, geo)

    return struc


def add_residue(
    structure: Structure, residue: Union[Geo, str], phi=-120, psi_im1=140, omega=-370
) -> Structure:
    """Add one residue to a peptide Structure using a residue code or a prepared geometry object."""

    if isinstance(residue, Geo):
        geo = residue
    elif isinstance(residue, str):
        geo = geometry(residue)
        geo.phi = phi
        geo.psi_im1 = psi_im1
        if omega > -361:
            geo.omega = omega
    else:
        raise ValueError("Invalid residue argument:", residue)

    return add_residue_from_geo(structure, geo)


def make_structure(
    AA_chain: str, phi: List[float], psi_im1: List[float], omega: Optional[List] = None
) -> Structure:
    """Build a peptide Structure from a sequence and residue-by-residue backbone angles."""
    geo = geometry(AA_chain[0])
    struc = initialize_res(geo)

    if omega is None or not len(omega):
        for i in range(1, len(AA_chain)):
            AA = AA_chain[i]
            add_residue(struc, AA, phi[i - 1], psi_im1[i - 1])
    else:
        for i in range(1, len(AA_chain)):
            AA = AA_chain[i]
            add_residue(struc, AA, phi[i - 1], psi_im1[i - 1], omega[i - 1])

    return struc


def make_structure_from_geos(geos: List[Geo]) -> Structure:
    """Build a peptide Structure from a list of prepared residue geometry objects."""
    model_structure = initialize_res(geos[0])
    for i in range(1, len(geos)):
        add_residue(model_structure, geos[i])

    return model_structure


def add_terminal_OXT(structure: Structure, C_OXT_length: float = 1.23) -> Structure:
    """Add the terminal OXT atom to the final residue of a completed peptide Structure."""

    rad = 180.0 / math.pi

    # obtain last residue infomation
    resRef = getReferenceResidue(structure)
    N_resRef = resRef["N"]
    CA_resRef = resRef["CA"]
    C_resRef = resRef["C"]
    O_resRef = resRef["O"]

    n_vec = N_resRef.get_vector()
    ca_vec = CA_resRef.get_vector()
    c_vec = C_resRef.get_vector()
    o_vec = O_resRef.get_vector()

    # geometry to bring together residue
    CA_C_OXT_angle = calc_angle(ca_vec, c_vec, o_vec) * rad
    N_CA_C_O_diangle = calc_dihedral(n_vec, ca_vec, c_vec, o_vec) * rad
    N_CA_C_OXT_diangle = N_CA_C_O_diangle - 180.0
    if N_CA_C_O_diangle < 0:
        N_CA_C_OXT_diangle = N_CA_C_O_diangle + 180.0

    # OXT atom creation
    OXT_coord = calculateCoordinates(
        N_resRef, CA_resRef, C_resRef, C_OXT_length, CA_C_OXT_angle, N_CA_C_OXT_diangle
    )
    OXT = Atom("OXT", OXT_coord, 0.0, 1.0, " ", "OXT", 0, "O")

    # modify last residue of the structure to contain the OXT atom
    resRef.add(OXT)
    return structure


class Geometry:
    """Notebook-local namespace compatible with `from PeptideBuilder import Geometry`."""
    Geo = Geo
    GlyGeo = GlyGeo
    AlaGeo = AlaGeo
    SerGeo = SerGeo
    CysGeo = CysGeo
    ValGeo = ValGeo
    IleGeo = IleGeo
    LeuGeo = LeuGeo
    ThrGeo = ThrGeo
    ArgGeo = ArgGeo
    LysGeo = LysGeo
    AspGeo = AspGeo
    AsnGeo = AsnGeo
    GluGeo = GluGeo
    GlnGeo = GlnGeo
    MetGeo = MetGeo
    HisGeo = HisGeo
    ProGeo = ProGeo
    PheGeo = PheGeo
    TyrGeo = TyrGeo
    TrpGeo = TrpGeo
    geometry = staticmethod(geometry)


class PeptideBuilder:
    """Notebook-local namespace compatible with `import PeptideBuilder`."""
    calculateCoordinates = staticmethod(calculateCoordinates)
    initialize_res = staticmethod(initialize_res)
    add_residue = staticmethod(add_residue)
    add_residue_from_geo = staticmethod(add_residue_from_geo)
    add_terminal_OXT = staticmethod(add_terminal_OXT)
    make_extended_structure = staticmethod(make_extended_structure)
    make_structure = staticmethod(make_structure)
    make_structure_from_geos = staticmethod(make_structure_from_geos)


# %% Cell rotamer-template-peptide-builder
AA1_TO_ROTAMER = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIE",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


def structure_to_simple_dataframe(structure):
    """Function:
        Convert the local peptide Structure into the pandas atom table used by this
        notebook.

    Parameters
    ----------
    structure : Structure
        Local peptide structure to convert.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table with one atom per row and the notebook atom columns.
    """
    rows = []
    anum = 1
    for model in structure.child_list:
        for chain in model.child_list:
            for residue in chain.child_list:
                resid = residue.get_id()[1]
                for atom in residue.child_list:
                    x, y, z = atom.coord
                    rows.append(
                        {
                            "atom": "ATOM",
                            "anum": anum,
                            "atom_name": atom.name,
                            "aa_name": residue.resname,
                            "resid": resid,
                            "x": float(x),
                            "y": float(y),
                            "z": float(z),
                        }
                    )
                    anum += 1
    return pd.DataFrame(
        rows, columns=["atom", "anum", "atom_name", "aa_name", "resid", "x", "y", "z"]
    )


def read_rotamer_template(
    residue_name: str, rotamer_dir: Union[str, Path] = DEFAULT_ROTAMER_DIR
) -> pd.DataFrame:
    """Function:
        Read one amino-acid template, including side-chain hydrogens, from
        RotamerLibrary.

    Parameters
    ----------
    residue_name : str
        Residue name used by RotamerLibrary, such as ALA or NALA.
    rotamer_dir : str or pathlib.Path
        Folder containing RotamerLibrary templates.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the rotamer-template atoms.
    """
    library_path = Path(rotamer_dir).expanduser()
    if not library_path.exists():
        raise FileNotFoundError(
            f"RotamerLibrary folder does not exist: {library_path.resolve()}. "
            "Check that the entire RotamerLibrary folder is available."
        )
    if not library_path.is_dir():
        raise NotADirectoryError(
            f"RotamerLibrary path is not a folder: {library_path.resolve()}"
        )

    template_file = library_path / residue_name
    if not template_file.is_file():
        raise FileNotFoundError(
            f"Rotamer file does not exist: {template_file.resolve()}. "
            "Check the residue files inside RotamerLibrary."
        )
    return _read_simple_pdb(str(template_file))


def _kabsch_transform(moving_points, target_points):
    """Function:
        Find the best rigid rotation and translation between two matched groups of 3D
        points.

    Parameters
    ----------
    moving_points : numpy.ndarray with shape (n, 3)
        Coordinates to move.
    target_points : numpy.ndarray with shape (n, 3)
        Matched target coordinates.

    Returns
    -------
    tuple
        Tuple (rotation, translation), numpy arrays with shapes (3, 3) and (3,).
    """
    moving_points = np.asarray(moving_points, dtype=np.float64)
    target_points = np.asarray(target_points, dtype=np.float64)
    moving_center = moving_points.mean(axis=0)
    target_center = target_points.mean(axis=0)

    moving0 = moving_points - moving_center
    target0 = target_points - target_center
    covariance = moving0.T @ target0
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = target_center - moving_center @ rotation.T
    return rotation, translation


def _apply_transform(coords, rotation, translation):
    """Function:
        Apply one rotation matrix and translation vector to a group of coordinates.

    Parameters
    ----------
    coords : numpy.ndarray with shape (n, 3)
        Coordinates to transform.
    rotation : numpy.ndarray with shape (3, 3)
        Rotation matrix.
    translation : numpy.ndarray with shape (3,)
        Translation vector.

    Returns
    -------
    numpy.ndarray with shape (n, 3)
        Transformed coordinates.
    """
    coords = np.asarray(coords, dtype=np.float64)
    return coords @ rotation.T + translation


def rebuild_terminal_carbonyl_for_beta(structure, geo, psi_im1):
    """Function:
        Move the terminal carbonyl oxygen into the requested beta-strand geometry.

    Parameters
    ----------
    structure : Structure
        Local peptide structure whose final carbonyl oxygen is moved.
    geo : Geo
        Geometry settings for the final residue.
    psi_im1 : float or list of floats in degrees
        Backbone psi angle or angles.

    Returns
    -------
    None
        Modifies the terminal O coordinate inside structure.
    """
    terminal_res = getReferenceResidue(structure)
    ghost_n_coord = calculateCoordinates(
        terminal_res["N"],
        terminal_res["CA"],
        terminal_res["C"],
        geo.peptide_bond,
        geo.CA_C_N_angle,
        psi_im1,
    )
    ghost_n = Atom("N", ghost_n_coord, 0.0, 0.0, " ", "N", 0, "N")
    terminal_res["O"].set_coord(
        calculateCoordinates(
            ghost_n,
            terminal_res["CA"],
            terminal_res["C"],
            geo.C_O_length,
            geo.CA_C_O_angle,
            180.0,
        )
    )


def build_target_beta_backbone(
    sequence: str, angles: Sequence[float]
) -> pd.DataFrame:
    """Function:
        Build the requested amino-acid sequence with beta-strand backbone angles and
        builder CB directions.

    Parameters
    ----------
    sequence : str
        Peptide sequence using one-letter amino-acid codes.
    angles : list or tuple [phi, psi] in degrees
        Beta-strand backbone angles.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the unhydrogenated peptide backbone and builder CB atoms.
    """
    geos = []
    for i, aa in enumerate(sequence):
        geo = Geometry.geometry(aa)
        if i > 0:
            geo.phi = angles[0]
            geo.psi_im1 = angles[1]
        geos.append(geo)
    structure = PeptideBuilder.make_structure_from_geos(geos)
    rebuild_terminal_carbonyl_for_beta(structure, geos[-1], angles[1])
    PeptideBuilder.add_terminal_OXT(structure)
    return structure_to_simple_dataframe(structure)


# %% Cell polyA-transplant-builder
# Clean peptide builder: poly-A/G backbone + rotamer side-chain transplant

BACKBONE_TEMPLATE_ATOMS = {
    "N",
    "H",
    "H1",
    "H2",
    "H3",
    "CA",
    "HA",
    "HA2",
    "HA3",
    "C",
    "O",
    "OXT",
}


def _atom_xyz_from_residue(residue_df, atom_name):
    """Function:
        Return the xyz coordinate of one named atom from one residue.

    Parameters
    ----------
    residue_df : pandas.DataFrame
        Atom rows for one residue only. The table must contain the columns `atom_name`,
        `x`, `y`, and `z`; other columns are allowed.
    atom_name : str
        Exact PDB atom name to search for, such as `CA`, `CB`, or `O`. The search is
        case-sensitive and does not change the table.

    Returns
    -------
    numpy.ndarray or None
        Coordinate [x, y, z] from the first matching row. Returns None when the atom is
        absent.
    """
    hit = residue_df[residue_df["atom_name"] == atom_name]
    if hit.empty:
        return None
    return hit.iloc[0][["x", "y", "z"]].to_numpy(dtype=np.float64)


def _append_atom(rows, atom_name, aa_name, resid, xyz):
    """Function:
        Append one atom record to a list that will later become a peptide DataFrame.

    Parameters
    ----------
    rows : list of dict
        Atom records that will become a pandas DataFrame.
    atom_name : str
        Name of the atom (e.g., 'CA', 'CB', 'N').
    aa_name : str
        Three- or four-character residue name.
    resid : int
        Residue number.
    xyz : array-like of float, shape (3,)
        Cartesian coordinates [x, y, z] of the atom.

    Returns
    -------
    None
        Appends one atom dictionary to rows.
    """
    rows.append(
        {
            "atom": "ATOM",
            "anum": 0,
            "atom_name": atom_name,
            "aa_name": aa_name,
            "resid": int(resid),
            "x": float(xyz[0]),
            "y": float(xyz[1]),
            "z": float(xyz[2]),
        }
    )


def Add_backbone_hydrogen(
    backbone_df: pd.DataFrame, sequence: str
) -> pd.DataFrame:
    # Add only peptide backbone and alpha hydrogens to an Ala/Gly scaffold.
    """Function:
        Add only backbone and alpha-carbon hydrogens to the peptide backbone.

    Parameters
    ----------
    backbone_df : pandas.DataFrame
        Peptide backbone atom coordinates.
    sequence : str
        Peptide sequence in one-letter amino-acid codes.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the backbone atoms and backbone hydrogens.
    """
    df = backbone_df.copy()
    sequence = sequence.upper().strip()
    additions = []
    residues = sorted(df["resid"].unique())
    residue_lookup = {resid: df[df["resid"] == resid].copy() for resid in residues}

    for i, resid in enumerate(residues):
        res = residue_lookup[resid]
        original_aa = sequence[i]
        aa_name = res.iloc[0]["aa_name"]
        n = _atom_xyz_from_residue(res, "N")
        ca = _atom_xyz_from_residue(res, "CA")
        c = _atom_xyz_from_residue(res, "C")
        cb = _atom_xyz_from_residue(res, "CB")

        if n is not None and ca is not None:
            if i == 0 and c is not None:
                if _atom_xyz_from_residue(res, "H") is None:
                    _append_atom(
                        additions,
                        "H",
                        aa_name,
                        resid,
                        place_fourth_atom(c, ca, n, 1.01, 109.5, 120.0),
                    )
                if _atom_xyz_from_residue(res, "H2") is None:
                    _append_atom(
                        additions,
                        "H2",
                        aa_name,
                        resid,
                        place_fourth_atom(c, ca, n, 1.01, 109.5, -120.0),
                    )
            elif original_aa != "P" and _atom_xyz_from_residue(res, "H") is None:
                prev = residue_lookup[residues[i - 1]]
                prev_c = _atom_xyz_from_residue(prev, "C")
                if prev_c is not None:
                    h = _place_opposite_bisector(n, [prev_c, ca], 1.01)
                    if h is not None:
                        _append_atom(additions, "H", aa_name, resid, h)

        if n is not None and ca is not None and c is not None:
            if original_aa == "G":
                if _atom_xyz_from_residue(res, "HA2") is None:
                    _append_atom(
                        additions,
                        "HA2",
                        aa_name,
                        resid,
                        place_fourth_atom(n, c, ca, 1.09, 109.5, 120.0),
                    )
                if _atom_xyz_from_residue(res, "HA3") is None:
                    _append_atom(
                        additions,
                        "HA3",
                        aa_name,
                        resid,
                        place_fourth_atom(n, c, ca, 1.09, 109.5, -120.0),
                    )
            elif _atom_xyz_from_residue(res, "HA") is None:
                neighbors = [n, c] if cb is None else [n, c, cb]
                ha = _place_opposite_bisector(ca, neighbors, 1.09)
                if ha is not None:
                    _append_atom(additions, "HA", aa_name, resid, ha)

    if additions:
        df = pd.concat([df, pd.DataFrame(additions)], ignore_index=True)

    df["_order"] = range(len(df))
    df["_added_h"] = (df["anum"] == 0).astype(int)
    df = df.sort_values(["resid", "_added_h", "_order"]).drop(
        columns=["_order", "_added_h"]
    )
    df = df.reset_index(drop=True)
    df["anum"] = range(1, len(df) + 1)
    return df


def _cb_local_frame(residue_df):
    """Function:
        Build a local coordinate frame around CA and CB for side-chain transplantation.

    Parameters
    ----------
    residue_df : pandas.DataFrame
        Atom rows belonging to one residue.

    Returns
    -------
    numpy.ndarray with shape (3, 3)
        Local CA-CB coordinate frame.
    """
    ca = _atom_xyz_from_residue(residue_df, "CA")
    cb = _atom_xyz_from_residue(residue_df, "CB")
    if ca is None or cb is None:
        raise ValueError("Need CA and CB atoms to build a side-chain transplant frame.")

    x_axis = unit_vector(cb - ca, "CA-CB vector")
    y_axis = None
    for anchor_name in ("N", "C"):
        anchor = _atom_xyz_from_residue(residue_df, anchor_name)
        if anchor is None:
            continue
        anchor_vec = anchor - ca
        anchor_vec = anchor_vec - np.dot(anchor_vec, x_axis) * x_axis
        if np.linalg.norm(anchor_vec) >= 1e-8:
            y_axis = unit_vector(anchor_vec, f"{anchor_name}-CA plane vector")
            break

    if y_axis is None:
        raise ValueError("Need N or C atom to orient the side-chain transplant frame.")

    z_axis = unit_vector(np.cross(x_axis, y_axis), "side-chain frame normal")
    y_axis = unit_vector(np.cross(z_axis, x_axis), "side-chain frame y-axis")
    return np.column_stack([x_axis, y_axis, z_axis])


def Transplant(
    backbone_df: pd.DataFrame, sequence: str,
    rotamer_dir: Union[str, Path] = DEFAULT_ROTAMER_DIR,
    residue_map: Optional[dict] = None
) -> pd.DataFrame:
    # Keep the builder backbone/CB atoms, then add side-chain atoms from RotamerLibrary.
    """Function:
        Transplant side-chain from RotamerLibrary based on backbone and CB coordinates.

    Parameters
    ----------
    backbone_df : pandas.DataFrame
        Peptide backbone atom coordinates.
    sequence : str
        Peptide sequence (one-letter amino acid codes).
    rotamer_dir : str or pathlib.Path
        Directory path containing RotamerLibrary template files.
    residue_map : dict or None, default=None
        Mapping from 1-letter amino acid codes to RotamerLibrary residue names.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the full peptide with transplanted side-chain atoms.
    """
    sequence = sequence.upper().strip()
    residue_map = dict(AA1_TO_ROTAMER if residue_map is None else residue_map)
    df = backbone_df.copy()
    sidechain_parts = []

    for resid, aa in enumerate(sequence, start=1):
        residue_name = residue_map[aa]
        df.loc[df["resid"] == resid, "aa_name"] = residue_name
        keep_atoms = BACKBONE_TEMPLATE_ATOMS | {"CB"}
        df = df[(df["resid"] != resid) | df["atom_name"].isin(keep_atoms)].copy()

        if aa == "G":
            continue

        target_res = df[df["resid"] == resid].copy()
        template = read_rotamer_template(residue_name, rotamer_dir=rotamer_dir).copy()

        template_ca = _atom_xyz_from_residue(template, "CA")
        template_cb = _atom_xyz_from_residue(template, "CB")
        target_ca = _atom_xyz_from_residue(target_res, "CA")
        target_cb = _atom_xyz_from_residue(target_res, "CB")
        if (
            template_ca is None
            or template_cb is None
            or target_ca is None
            or target_cb is None
        ):
            raise ValueError(
                f"Need CA and CB atoms to transplant residue {residue_name} "
                f"at resid {resid}."
            )

        template_cb_length = np.linalg.norm(template_cb - template_ca)
        target_cb_direction = unit_vector(target_cb - target_ca, "builder CA-CB vector")
        placed_cb = target_ca + target_cb_direction * template_cb_length
        cb_mask = (df["resid"] == resid) & (df["atom_name"] == "CB")
        df.loc[cb_mask, ["x", "y", "z"]] = placed_cb
        target_res.loc[target_res["atom_name"] == "CB", ["x", "y", "z"]] = placed_cb

        rotation = _cb_local_frame(target_res) @ _cb_local_frame(template).T

        side_atoms = template[
            ~template["atom_name"].isin(BACKBONE_TEMPLATE_ATOMS | {"CB"})
        ].copy()
        if side_atoms.empty:
            continue

        coords = side_atoms[["x", "y", "z"]].to_numpy(dtype=np.float64)
        side_atoms[["x", "y", "z"]] = (coords - template_cb) @ rotation.T + placed_cb
        side_atoms["resid"] = resid
        side_atoms["aa_name"] = residue_name
        sidechain_parts.append(side_atoms)

    if sidechain_parts:
        df = pd.concat([df] + sidechain_parts, ignore_index=True)

    df["_order"] = range(len(df))
    df = (
        df.sort_values(["resid", "_order"])
        .drop(columns=["_order"])
        .reset_index(drop=True)
    )
    df["anum"] = range(1, len(df) + 1)
    return df


def build_single_peptide_with_local_peptidebuilder(
    sequence: str, angles: Sequence[float],
    rotamer_dir: Union[str, Path] = DEFAULT_ROTAMER_DIR,
    residue_map: Optional[dict] = None
) -> pd.DataFrame:
    """Function:
        Build one peptide DataFrame without writing an intermediate PDB file.

    Parameters
    ----------
    sequence : str
        Peptide sequence using one-letter amino-acid codes.
    angles : list or tuple [phi, psi] in degrees
        Backbone angles used for the full strand.
    rotamer_dir : str or pathlib.Path
        Folder containing RotamerLibrary templates.
    residue_map : dict or None
        Maps one-letter amino-acid codes to RotamerLibrary residue names.

    Returns
    -------
    pandas.DataFrame
        Generated peptide atoms with backbone hydrogens and side chains.
    """
    sequence = sequence.upper().strip()

    scaffold_sequence = sequence

    backbone_df = build_target_beta_backbone(scaffold_sequence, angles)

    backbone_df = Add_backbone_hydrogen(backbone_df, sequence)
    peptide_df = Transplant(
        backbone_df, sequence, rotamer_dir=rotamer_dir, residue_map=residue_map
    )
    return peptide_df.reset_index(drop=True)


# %% Cell 990c2583-b6fa-4d1d-ab63-c54ad39ec366
def balance_oxygen_z_roll(df):
    """Function:
        Roll one strand around y so the selected backbone oxygen atoms are balanced
        around the x-y plane.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table after a y-axis roll that balances selected O atoms around the
        x-y plane.
    """
    peptide_residues = (
        df[
            (df["atom_name"] == "CA")
            & (df["aa_name"] != "ACE")
            & (df["aa_name"] != "NME")
        ]["resid"]
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    if len(peptide_residues) < 3:
        return df
    if len(peptide_residues) % 2 == 1:
        peptide_residues = peptide_residues[:-1]
    residue_index = {resid: i + 1 for i, resid in enumerate(peptide_residues)}
    selected_oxygen = df[
        (df["atom_name"] == "O") & (df["resid"].isin(peptide_residues))
    ].copy()
    selected_oxygen["_peptide_index"] = selected_oxygen["resid"].map(residue_index)
    selected_oxygen = selected_oxygen.sort_values("_peptide_index")
    if len(selected_oxygen) < 2:
        return df

    x = selected_oxygen["x"].to_numpy(dtype=np.float64)
    z = selected_oxygen["z"].to_numpy(dtype=np.float64)
    x0 = x - x.mean()
    z0 = z - z.mean()
    s_xx = float(np.dot(x0, x0))
    s_zz = float(np.dot(z0, z0))
    s_xz = float(np.dot(x0, z0))
    if s_xx + s_zz < 1e-10:
        return df

    theta0 = 0.5 * np.arctan2(2.0 * s_xz, s_xx - s_zz)
    candidates = [theta0, theta0 + np.pi / 2.0, theta0 - np.pi / 2.0]

    def z_spread_after_roll(theta):
        """Function:
            Calculate the oxygen z-coordinate variance after one trial roll.

        Parameters
        ----------
        theta : float
            Trial rotation angle around the y axis, in radians.

        Returns
        -------
        float
            Variance of the rotated oxygen z coordinates.
        """
        return float(np.var(-np.sin(theta) * x0 + np.cos(theta) * z0))

    theta = min(candidates, key=z_spread_after_roll)
    theta = min([theta, theta + np.pi, theta - np.pi], key=abs)
    df_roll = rotation_coordinates(df, "y", theta)

    gly_ha2 = df_roll["aa_name"].astype(str).str.endswith("GLY") & (
        df_roll["atom_name"] == "HA2"
    )
    center_atoms = df_roll[
        ((df_roll["atom_name"].isin(["N", "CA", "C", "CB"])) | gly_ha2)
        & (df_roll["aa_name"] != "ACE")
        & (df_roll["aa_name"] != "NME")
    ]
    if not center_atoms.empty:
        x_mid = 0.5 * (center_atoms["x"].min() + center_atoms["x"].max())
        df_roll["x"] = df_roll["x"] - x_mid
    return df_roll


def _center_peptide_x_after_roll(df):
    """Function:
        Move a rolled peptide back to the center of the x direction.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table translated so backbone, CB, and Gly HA2 atoms are centered in
        x.
    """
    df = df.copy()
    gly_ha2 = df["aa_name"].astype(str).str.endswith("GLY") & (df["atom_name"] == "HA2")
    center_atoms = df[
        ((df["atom_name"].isin(["N", "CA", "C", "CB"])) | gly_ha2)
        & (df["aa_name"] != "ACE")
        & (df["aa_name"] != "NME")
    ]
    if not center_atoms.empty:
        x_mid = 0.5 * (center_atoms["x"].min() + center_atoms["x"].max())
        df["x"] = df["x"] - x_mid
    return df


def _middle_cb_residue_ids(df):
    """Function:
        Select the central residue or residues used for the parallel-strand CA-CB angle
        check.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.

    Returns
    -------
    list of int residue numbers
        One center residue for odd length or two for even length.
    """
    peptide_residues = (
        df[
            (df["atom_name"] == "CA")
            & (df["aa_name"] != "ACE")
            & (df["aa_name"] != "NME")
        ]["resid"]
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    cb_residues = set(df[df["atom_name"] == "CB"]["resid"].tolist())
    if not peptide_residues or not cb_residues:
        return []

    center_index = 0.5 * (len(peptide_residues) - 1)
    candidates = [
        (abs(index - center_index), resid)
        for index, resid in enumerate(peptide_residues)
        if resid in cb_residues
    ]
    if not candidates:
        return []
    nearest_distance = min(distance for distance, _ in candidates)
    return [
        resid
        for distance, resid in candidates
        if abs(distance - nearest_distance) < 1e-8
    ]


def ca_cb_tilt_angles(df, residue_ids=None):
    """Function:
        Measure selected CA-CB bond angles from the nearest positive or negative z
        direction.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    residue_ids : list of int or None
        Residues to measure. None uses every residue that has a CB atom.

    Returns
    -------
    numpy.ndarray with shape (m,)
        CA-CB tilt angles in degrees from the nearest +z/-z direction.
    """
    angles = []
    if residue_ids is None:
        residue_ids = (
            df[(df["aa_name"] != "ACE") & (df["aa_name"] != "NME")]["resid"]
            .drop_duplicates()
            .tolist()
        )
    for resid in residue_ids:
        residue = df[df["resid"] == resid]
        ca = residue[residue["atom_name"] == "CA"]
        cb = residue[residue["atom_name"] == "CB"]
        if ca.empty or cb.empty:
            continue
        ca_xyz = ca[["x", "y", "z"]].iloc[0].to_numpy(dtype=np.float64)
        cb_xyz = cb[["x", "y", "z"]].iloc[0].to_numpy(dtype=np.float64)
        v = cb_xyz - ca_xyz
        v_norm = np.linalg.norm(v)
        if v_norm < 1e-8:
            continue
        cos_to_z = abs(v[2] / v_norm)
        angles.append(np.degrees(np.arccos(np.clip(cos_to_z, -1.0, 1.0))))
    return np.asarray(angles, dtype=np.float64)


def restrict_parallel_ca_cb_tilt(df, max_cb_tilt_deg=10.0):
    """Function:
        Apply the smallest y rotation needed to keep only the middle CA-CB bond or bonds
        within the angle limit.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    max_cb_tilt_deg : float or None
        Maximum allowed middle CA-CB tilt in degrees.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table after the smallest y-axis roll that satisfies the
        center-angle limit.
    """
    middle_residue_ids = _middle_cb_residue_ids(df)
    angles = ca_cb_tilt_angles(df, residue_ids=middle_residue_ids)
    if len(angles) == 0 or float(angles.max()) <= max_cb_tilt_deg:
        return df

    best_df = df
    best_max_tilt = float(angles.max())
    feasible = []

    # Use a small y-axis roll only. This keeps the strand axis and packing direction
    # unchanged.
    for theta in np.linspace(-np.pi / 2.0, np.pi / 2.0, 721):
        rolled = rotation_coordinates(df, "y", float(theta))
        rolled = _center_peptide_x_after_roll(rolled)
        rolled_angles = ca_cb_tilt_angles(rolled, residue_ids=middle_residue_ids)
        if len(rolled_angles) == 0:
            continue
        max_tilt = float(rolled_angles.max())
        if max_tilt < best_max_tilt:
            best_max_tilt = max_tilt
            best_df = rolled
        if max_tilt <= max_cb_tilt_deg + 1e-6:
            feasible.append((abs(float(theta)), max_tilt, rolled))

    if feasible:
        feasible.sort(key=lambda item: (item[0], item[1]))
        return feasible[0][2]
    warnings.warn(
        f"Cannot make middle CA-CB angle within {max_cb_tilt_deg:.1f} "
        f"degree; best angle is {best_max_tilt:.2f} degree."
    )
    return best_df


def alignment(
    length: int, chains: int, peptide_df: pd.DataFrame,
    o_filename: Union[str, Path], order: int = 1,
    max_cb_tilt_deg: Optional[float] = None
) -> pd.DataFrame:
    """Function:
        Align a peptide DataFrame and write one aligned checking PDB file.

    Parameters
    ----------
    length : int
        Number of peptide residues, including caps when present.
    chains : int
        Number of peptide chains in the DataFrame, normally 1.
    peptide_df : pandas.DataFrame
        Peptide atoms after terminal processing.
    o_filename : str or pathlib.Path
        Output filename for the aligned checking PDB.
    order : int
        Selects the initial strand-axis direction during alignment.
    max_cb_tilt_deg : float or None
        Maximum allowed middle CA-CB tilt in degrees.

    Returns
    -------
    pandas.DataFrame
        Aligned peptide atoms. The function also writes the checking PDB.
    """
    df = peptide_df.copy().reset_index(drop=True)

    if df[df["resid"] == 1]["aa_name"].iloc[0] != "ACE":
        residue_window = df["resid"] <= 10
    else:
        residue_window = (
            (df["resid"] > 1) & (df["resid"] <= 11) & (df["aa_name"] != "NME")
        )
    gly_ha2 = df["aa_name"].astype(str).str.endswith("GLY") & (df["atom_name"] == "HA2")
    filtered_df = df[
        residue_window & (df["atom_name"].isin(["C", "CA", "CB"]) | gly_ha2)
    ]

    # filtered_df = df[(df['resid'] <= length * chains) & (df['atom_name'].isin(['O',
    # 'C']))]
    data = filtered_df[["x", "y", "z"]].to_numpy()

    # regular grid covering the domain of the data
    X_min = filtered_df["x"].min()
    X_max = filtered_df["x"].max()
    Y_min = filtered_df["y"].min()
    Y_max = filtered_df["y"].max()
    X, Y = np.meshgrid(np.arange(X_min, X_max, 0.6), np.arange(Y_min, Y_max, 0.5))
    XX = X.flatten()
    YY = Y.flatten()

    # ################# fitting sheet plane and get the plane normal vector
    # #################
    # order = 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        # data[:,0] = x-column, data[:,1] = y-column, data[:,2] = z-column
        A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
        C, _, _, _ = np.linalg.lstsq(A, data[:, 2], rcond=None)  # coefficients

        # evaluate it on grid
        Z = C[0] * X + C[1] * Y + C[2]

        # or expressed using matrix/vector product
        # Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[
            np.ones(data.shape[0]),
            data[:, :2],
            np.prod(data[:, :2], axis=1),
            data[:, :2] ** 2,
        ]
        C, _, _, _ = np.linalg.lstsq(A, data[:, 2], rcond=None)
        print(A)
        # evaluate it on a grid
        Z = np.dot(
            np.c_[np.ones(XX.shape), XX, YY, XX * YY, XX**2, YY**2], C
        ).reshape(X.shape)

    # ###########calculate the angles between z-axis and the plane normal vector
    # ###################
    nz = [0, 0, 1]
    # sheets plane normal vector = (C[0], C[1], -1)
    # rotation on x-y plane (rotate by z-axis). v1 = (0,0), *** no need to rotate by
    # z-axis ***
    # v1=[C[0],C[1]]
    # v2=[nz[0],nz[1]]
    # tz = rotation_angle(v1,v2) # rotation from v1 to v2

    # rotation on x-z plane (rotate by y-axis)
    v1 = [C[0], -1]
    v2 = [nz[0], nz[2]]
    ty = rotation_angle(v1, v2)  # rotation from v1 to v2

    # rotation on y-z plane (rotate by x-axis)
    v1 = [C[1], -1]
    v2 = [nz[1], nz[2]]
    tx = rotation_angle(v1, v2)  # rotation from v1 to v2

    # ###########calculate rotation matrix and rotate the plane normal to z-axis
    # ###################
    Rx_M = np.array(
        [[1, 0, 0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]]
    )
    Ry_M = np.array(
        [[np.cos(ty), 0, np.sin(ty)], [0, 1, 0], [-np.sin(ty), 0, np.cos(ty)]]
    )
    # Rz_M = np.array([[np.cos(ty), 0, np.sin(ty)], [0, 1, 0], [-np.sin(ty), 0,
    # np.cos(ty)]])

    data = df[["x", "y", "z"]].to_numpy()
    df_temp = df.copy()

    for i, xyz in enumerate(data):
        xyz = Rx_M @ xyz
        xyz = Ry_M @ xyz

        df_temp.at[i, "x"] = xyz[0]
        df_temp.at[i, "y"] = xyz[1]
        df_temp.at[i, "z"] = xyz[2]

    # ############################ fit each peptide parallel to y-axis and move the
    # sheets to center ###############
    vectors = []

    # filtered_df = df_temp[(df_temp['resid']>i*length)&(df_temp['resid']<=(i+1)*length)
    # &(df_temp['atom_name'].isin(['C','CA','N']))]
    filtered_df = df_temp[df_temp["atom_name"].isin(["CA"])]
    points = filtered_df[["x", "y", "z"]].to_numpy()

    #       k1, b1, k2, b2 = linear_fitting_3D_points(points)

    #         p = 1
    #         m = k1 * p
    #         n = k2 * p
    m, n, p = linear_fitting(points)
    vector = [m, n, p]  # sheet direction vector
    vectors.append(vector)

    vectors = np.array(vectors)
    vector = np.mean(vectors, axis=0)

    # ############################ align peptide to y-axis and move center to origin
    # ###############
    ny = [0, 1, 0]
    v1 = [vector[0], vector[1]]
    v2 = [ny[0], ny[1]]
    tz = rotation_angle(v1, v2)  # rotation from v1 to v2

    Rz_M = np.array(
        [[np.cos(tz), -np.sin(tz), 0], [np.sin(tz), np.cos(tz), 0], [0, 0, 1]]
    )

    data = df_temp[["x", "y", "z"]].to_numpy()

    for i, xyz in enumerate(data):
        xyz = Rz_M @ xyz
        df_temp.at[i, "x"] = xyz[0]
        df_temp.at[i, "y"] = xyz[1]
        df_temp.at[i, "z"] = xyz[2]

    x_sum = 0
    y_sum = 0
    z_sum = 0
    data_for_center = df_temp[
        (df_temp["atom_name"].isin(["C", "CA", "N"]))
        & (df_temp["aa_name"] != "ACE")
        & (df_temp["aa_name"] != "NME")
    ][["x", "y", "z"]].to_numpy()
    for i, xyz in enumerate(data_for_center):
        x_sum += xyz[0]
        y_sum += xyz[1]
        z_sum += xyz[2]

    x_sum /= i + 1
    y_sum /= i + 1
    z_sum /= i + 1
    df_temp["x"] = df_temp["x"] - x_sum
    df_temp["y"] = df_temp["y"] - y_sum
    df_temp["z"] = df_temp["z"] - z_sum

    if df[df["resid"] == 1]["aa_name"].iloc[0] != "ACE":
        res_1_CA_x = df_temp[(df_temp["resid"] == 1) & (df_temp["atom_name"] == "CA")][
            "x"
        ].iloc[0]
        res_2_CA_x = df_temp[(df_temp["resid"] == 2) & (df_temp["atom_name"] == "CA")][
            "x"
        ].iloc[0]
        res_1_CA_y = df_temp[(df_temp["resid"] == 1) & (df_temp["atom_name"] == "CA")][
            "y"
        ].iloc[0]
        res_2_CA_y = df_temp[(df_temp["resid"] == 2) & (df_temp["atom_name"] == "CA")][
            "y"
        ].iloc[0]
    else:
        res_1_CA_x = df_temp[(df_temp["resid"] == 2) & (df_temp["atom_name"] == "CA")][
            "x"
        ].iloc[0]
        res_2_CA_x = df_temp[(df_temp["resid"] == 3) & (df_temp["atom_name"] == "CA")][
            "x"
        ].iloc[0]
        res_1_CA_y = df_temp[(df_temp["resid"] == 2) & (df_temp["atom_name"] == "CA")][
            "y"
        ].iloc[0]
        res_2_CA_y = df_temp[(df_temp["resid"] == 3) & (df_temp["atom_name"] == "CA")][
            "y"
        ].iloc[0]

    if res_1_CA_x > res_2_CA_x:
        ty = np.pi / 2
    else:
        ty = -np.pi / 2

    df_temp = rotation_coordinates(df_temp, "y", ty)

    if res_1_CA_y > res_2_CA_y:
        tz = np.pi
        df_temp = rotation_coordinates(df_temp, "z", tz)

    # The beta strand is twisted, not a perfect plane. Center the twist around
    # the y-z midplane so roughly half of the strand envelope is at x < 0 and
    # half is at x > 0 when viewed along the strand axis.
    gly_ha2 = df_temp["aa_name"].astype(str).str.endswith("GLY") & (
        df_temp["atom_name"] == "HA2"
    )
    spiral_center_atoms = df_temp[
        ((df_temp["atom_name"].isin(["N", "CA", "C", "CB"])) | gly_ha2)
        & (df_temp["aa_name"] != "ACE")
        & (df_temp["aa_name"] != "NME")
    ]
    if not spiral_center_atoms.empty:
        x_mid = 0.5 * (spiral_center_atoms["x"].min() + spiral_center_atoms["x"].max())
        df_temp["x"] = df_temp["x"] - x_mid

    df_temp = balance_oxygen_z_roll(df_temp)
    if max_cb_tilt_deg is not None:
        df_temp = restrict_parallel_ca_cb_tilt(df_temp, max_cb_tilt_deg=max_cb_tilt_deg)

    #     df_temp = df_temp[(df_temp['resid'] <= length)]

    ############# write out the transformed structure ###################
    with open(o_filename, "w") as f:
        for index, row in df_temp.iterrows():
            if len(row["aa_name"]) == 4:
                output = "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row["anum"],
                    row["atom_name"],
                    row["aa_name"],
                    row["resid"],
                    row["x"],
                    row["y"],
                    row["z"],
                )
            else:
                output = "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row["anum"],
                    row["atom_name"],
                    row["aa_name"],
                    row["resid"],
                    row["x"],
                    row["y"],
                    row["z"],
                )

            f.write(output)

    return df_temp


# %% Cell 54ce3e23-782f-45a3-8fc4-31e20d984620
def _packing_core_atoms(df):
    """Function:
        Select backbone, CB, and Gly HA2 atoms used to calculate packing centers.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.

    Returns
    -------
    pandas.DataFrame
        DataFrame subset containing N, CA, C, CB, and Gly HA2 atoms.
    """
    gly_ha2 = df["aa_name"].astype(str).str.endswith("GLY") & (df["atom_name"] == "HA2")
    core = df[
        ((df["atom_name"].isin(["N", "CA", "C", "CB"])) | gly_ha2)
        & (df["aa_name"] != "ACE")
        & (df["aa_name"] != "NME")
    ]
    return df if core.empty else core


def _packing_center(df):
    """Function:
        Calculate the mean xyz center of the atoms used for packing.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Mean [x, y, z] packing center.
    """
    return (
        _packing_core_atoms(df)[["x", "y", "z"]].to_numpy(dtype=np.float64).mean(axis=0)
    )


def _packing_axis_midpoint(df, coord):
    """Function:
        Calculate the midpoint between the minimum and maximum value on one coordinate
        axis.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    coord : str
        Coordinate column name 'x', 'y', or 'z'.

    Returns
    -------
    float
        Midpoint between the minimum and maximum coordinate.
    """
    values = _packing_core_atoms(df)[coord].to_numpy(dtype=np.float64)
    return 0.5 * (values.min() + values.max())


def _chain_center(df, chain_length, chain_index=0):
    """Function:
        Calculate the packing center of one selected chain in a packed sheet.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    chain_index : int
        Zero-based index of the chain to select.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Center of the selected chain.
    """
    return _packing_center(_chain_slice(df, chain_length, chain_index=chain_index))


def _translate_peptide(df, dx=0.0, dy=0.0, dz=0.0):
    """Function:
        Move every atom in a peptide by defined x, y, and z distances.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    dx : float
        Translation distance along x in Angstrom.
    dy : float
        Translation distance along y in Angstrom.
    dz : float
        Translation distance along z in Angstrom.

    Returns
    -------
    pandas.DataFrame
        Translated copy of df.
    """
    df = df.copy()
    df["x"] = df["x"] + dx
    df["y"] = df["y"] + dy
    df["z"] = df["z"] + dz
    return df


def _ca_trace(df):
    """Function:
        Select and order the CA atoms that define one peptide strand trace.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.

    Returns
    -------
    pandas.DataFrame
        Atom table containing CA rows sorted by residue number.
    """
    return df[
        (df["atom_name"] == "CA") & (df["aa_name"] != "ACE") & (df["aa_name"] != "NME")
    ].sort_values("resid")


def _chain_slice(df, chain_length, chain_index=0):
    """Function:
        Select one chain from a DataFrame containing several packed chains.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    chain_index : int
        Zero-based index of the chain to select.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the selected chain.
    """
    min_resid = int(df["resid"].min())
    start = min_resid + chain_index * chain_length
    stop = start + chain_length
    return df[(df["resid"] >= start) & (df["resid"] < stop)]


def _strand_axis(df, chain_length=None, chain_index=0):
    """Function:
        Fit the N-to-C direction of one strand from its ordered CA atoms.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    chain_length : int or None
        Use this value to select one chain when df contains packed chains.
    chain_index : int
        Zero-based index of the chain to select.

    Returns
    -------
    numpy.ndarray with shape (3,)
        Unit N-to-C strand direction.
    """
    if chain_length is not None:
        df = _chain_slice(df, chain_length, chain_index=chain_index)
    ca = _ca_trace(df)
    points = ca[["x", "y", "z"]].to_numpy(dtype=np.float64)
    if len(points) < 2:
        raise ValueError("Need at least two CA atoms to define a strand axis.")
    axis = linear_fitting(points)
    first_to_last = points[-1] - points[0]
    if np.dot(axis, first_to_last) < 0:
        axis = -axis
    return unit_vector(axis, "strand axis")


def _xy_angle_between(v_from, v_to):
    """Function:
        Calculate the signed rotation angle between two vectors after projection onto
        the x-y plane.

    Parameters
    ----------
    v_from : array-like vector
        Starting direction for the angle calculation.
    v_to : array-like vector
        Target direction for the angle calculation.

    Returns
    -------
    float
        Signed angle in radians in the x-y plane.
    """
    a = unit_vector(
        np.asarray([v_from[0], v_from[1]], dtype=np.float64), "moving xy strand axis"
    )
    b = unit_vector(
        np.asarray([v_to[0], v_to[1]], dtype=np.float64), "target xy strand axis"
    )
    return np.arctan2(cross2d(a, b), np.dot(a, b))


def _rotate_z_about_center(df, angle, center=None):
    """Function:
        Rotate a peptide around a z-axis passing through a defined center.

    Parameters
    ----------
    df : pandas.DataFrame
        Peptide atom table containing x, y, z and residue information.
    angle : float
        Rotation angle in radians.
    center : numpy.ndarray with shape (3,) or None
        Rotation center. None uses the calculated packing center.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table after rotation around z through center.
    """
    df = df.copy()
    coords = df[["x", "y", "z"]].to_numpy(dtype=np.float64)
    if center is None:
        center = _packing_center(df)
    center = np.asarray(center, dtype=np.float64)
    rot = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    df[["x", "y", "z"]] = (coords - center) @ rot.T + center
    return df


def _align_strand_xy(moving_df, target_df, same_direction=None, chain_length=None,
                     moving_chain_index=0, target_chain_index=0):
    """Function:
        Rotate a moving strand so its projected strand direction matches a target
        strand.

    Parameters
    ----------
    moving_df : pandas.DataFrame
        Strand or sheet coordinates that will be moved.
    target_df : pandas.DataFrame
        Fixed reference strand or sheet coordinates.
    same_direction : bool or None
        True aligns the same N-to-C direction, False aligns the opposite direction, and
        None chooses automatically.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    moving_chain_index : int
        Zero-based moving-chain index when moving_df contains several chains.
    target_chain_index : int
        Zero-based target-chain index when target_df contains several chains.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the rotated moving strand or sheet.
    """
    moving_axis = _strand_axis(
        moving_df, chain_length=chain_length, chain_index=moving_chain_index
    )
    target_axis = _strand_axis(
        target_df, chain_length=chain_length, chain_index=target_chain_index
    )
    if same_direction is None:
        desired_axis = (
            target_axis
            if np.dot(moving_axis[:2], target_axis[:2]) >= 0
            else -target_axis
        )
    else:
        desired_axis = target_axis if same_direction else -target_axis
    angle = _xy_angle_between(moving_axis, desired_axis)
    return _rotate_z_about_center(moving_df, angle)


def _register_antiparallel_ca_y(moving_df, target_df, residue_offset=0,
                                direction="minus"):
    """Function:
        Translate an antiparallel strand along y to match the requested CA registry.

    Parameters
    ----------
    moving_df : pandas.DataFrame
        Strand or sheet coordinates that will be moved.
    target_df : pandas.DataFrame
        Fixed reference strand or sheet coordinates.
    residue_offset : int
        Number of residues used for the registry offset.
    direction : str
        'minus' or 'plus' selects the registry-shift direction.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table translated to the antiparallel CA registry.
    """
    moving_ca = _ca_trace(moving_df)
    target_ca = _ca_trace(target_df)
    n = min(len(moving_ca), len(target_ca))
    if n == 0:
        return moving_df
    moving_y = moving_ca["y"].to_numpy(dtype=np.float64)[:n][::-1]
    target_y = target_ca["y"].to_numpy(dtype=np.float64)[:n]
    use_anchor_pair = False
    if residue_offset:
        residue_offset = min(abs(int(residue_offset)), n - 1)
        if direction == "plus":
            moving_y = moving_y[: n - residue_offset]
            target_y = target_y[residue_offset:n]
        else:
            moving_y = moving_y[residue_offset:]
            target_y = target_y[: len(moving_y)]
        use_anchor_pair = True
    dy = (
        float(target_y[0] - moving_y[0])
        if use_anchor_pair
        else float(np.mean(target_y - moving_y))
    )
    return _translate_peptide(moving_df, dy=dy)


def _register_parallel_ca_yz(moving_df, target_df, residue_offset=0, direction="minus"):
    """Function:
        Rigidly align a parallel strand to the target CA positions in y and z.

    Parameters
    ----------
    moving_df : pandas.DataFrame
        Strand or sheet coordinates that will be moved.
    target_df : pandas.DataFrame
        Fixed reference strand or sheet coordinates.
    residue_offset : int
        Number of residues used for the registry offset.
    direction : str
        'minus' or 'plus' selects the registry-shift direction.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table transformed to the parallel CA registry.
    """
    moving_ca = _ca_trace(moving_df)
    target_ca = _ca_trace(target_df)
    n = min(len(moving_ca), len(target_ca))
    if n == 0:
        return moving_df
    moving_points = moving_ca[["x", "y", "z"]].to_numpy(dtype=np.float64)[:n]
    target_points = target_ca[["x", "y", "z"]].to_numpy(dtype=np.float64)[:n]
    if residue_offset:
        residue_offset = min(abs(int(residue_offset)), n - 1)
        if direction == "plus":
            moving_points = moving_points[residue_offset:n]
            target_points = target_points[: len(moving_points)]
        else:
            moving_points = moving_points[: n - residue_offset]
            target_points = target_points[residue_offset:n]
    if len(moving_points) < 2:
        delta = target_points.mean(axis=0) - moving_points.mean(axis=0)
        return _translate_peptide(
            moving_df, dx=float(delta[0]), dy=float(delta[1]), dz=float(delta[2])
        )
    rotation, translation = _kabsch_transform(moving_points, target_points)
    registered = moving_df.copy()
    coords = registered[["x", "y", "z"]].to_numpy(dtype=np.float64)
    registered[["x", "y", "z"]] = _apply_transform(coords, rotation, translation)
    return registered


def _estimate_residue_spacing(single_peptite):
    """Function:
        Estimate the neighboring-residue spacing from the median CA-to-CA y distance.

    Parameters
    ----------
    single_peptite : pandas.DataFrame
        One aligned peptide strand used as a packing template.

    Returns
    -------
    float
        Median neighboring CA y-spacing in Angstrom.
    """
    ca = single_peptite[
        (single_peptite["atom_name"] == "CA")
        & (single_peptite["aa_name"] != "ACE")
        & (single_peptite["aa_name"] != "NME")
    ].sort_values("resid")
    y_coords = ca["y"].to_numpy(dtype=np.float64)
    if len(y_coords) < 2:
        return 3.3
    spacings = np.abs(np.diff(y_coords))
    spacings = spacings[spacings > 1e-6]
    if len(spacings) == 0:
        return 3.3
    return float(np.median(spacings))


def _write_packed_peptides(peptides, chain_length, format_flag, o_filename):
    """Function:
        Write a packed one-component peptide DataFrame to a PDB file.

    Parameters
    ----------
    peptides : pandas.DataFrame
        All packed peptide atom rows.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    format_flag : int
        0 selects PepAD format; 1 selects AMBER format.
    o_filename : str or pathlib.Path
        Output PDB file.

    Returns
    -------
    None
        Writes the packed PDB file.
    """
    with open(o_filename, "w") as f:
        last_resid = None
        for _, row in peptides.iterrows():
            resid = int(row["resid"])
            aa_name = (
                row["aa_name"][1:]
                if format_flag == 1 and len(row["aa_name"]) == 4
                else row["aa_name"]
            )
            if format_flag == 1 and last_resid is not None:
                if (last_resid % chain_length == 0) and (resid % chain_length == 1):
                    f.write("TER\n")
            output = "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                int(row["anum"]),
                row["atom_name"],
                aa_name,
                resid,
                float(row["x"]),
                float(row["y"]),
                float(row["z"]),
            )
            f.write(output)
            last_resid = resid


def packpep_geometric(single_peptite, chain_length, classes, shifts, chains_per_sheet,
                      format_flag, o_filename, core_residues="e", strand_dist=4.8,
                      sheet_dist=11.5, sheet_shift=0.0, registry_direction="minus"):
    """Function:
        Construct and pack two one-component sheets for steric-zipper Classes 1-10.

    Parameters
    ----------
    single_peptite : pandas.DataFrame
        One aligned peptide strand used as a packing template.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    classes : int
        Steric-zipper class number from 1 to 10.
    shifts : int or float
        Move sheet 2 along y by this number of residue spacings.
    chains_per_sheet : int
        Number of peptide strands placed in each sheet.
    format_flag : int
        0 selects PepAD format; 1 selects AMBER format.
    o_filename : str or pathlib.Path
        Output PDB file.
    core_residues : str
        'e' packs even-numbered residues inside and 'o' packs odd-numbered residues
        inside.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    sheet_dist : float
        Sheet-to-sheet distance along z in Angstrom.
    sheet_shift : float
        Move sheet 2 along x in units of half the strand distance.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    pandas.DataFrame containing both packed sheets
        Also writes o_filename.
    """
    class_axis = {
        1: "x",
        2: "none",
        3: "y",
        4: "z",
        5: "x",
        6: "none",
        7: "none",
        8: "z",
        9: "none",
        10: "x",
    }
    if classes not in class_axis:
        raise ValueError("classes must be an integer from 1 to 10")

    anum_max = int(single_peptite["anum"].max())
    resid_max = int(single_peptite["resid"].max())
    if registry_direction not in ["minus", "plus"]:
        raise ValueError("registry_direction must be 'minus' or 'plus'")
    residue_spacing = _estimate_residue_spacing(single_peptite)
    sheet_y_shift = float(shifts) * residue_spacing
    sheet_x_shift = float(sheet_shift) * float(strand_dist) * 0.5

    if classes in [5, 6, 7, 8]:
        anti_axis = "z" if classes in [5, 6] else "x"
        anti_pep = rotation_coordinates(single_peptite, anti_axis, np.pi)
        anti_pep = _align_strand_xy(anti_pep, single_peptite, same_direction=False)
        if classes in [5, 6]:
            register_offset = 1 if chain_length % 2 == 0 else 0
        elif classes in [7, 8]:
            register_offset = 1 if chain_length % 2 == 1 else 0
        else:
            register_offset = 0
        anti_pep = _register_antiparallel_ca_y(
            anti_pep,
            single_peptite,
            residue_offset=register_offset,
            direction=registry_direction,
        )
        center_delta = _packing_center(single_peptite) - _packing_center(anti_pep)
        anti_pep = _translate_peptide(anti_pep, center_delta[0], 0.0, center_delta[2])
    else:
        anti_pep = None

    if classes in [9, 10]:
        face_pep = rotation_coordinates(single_peptite, "y", np.pi)
        face_pep = _align_strand_xy(face_pep, single_peptite, same_direction=True)
        face_delta = _packing_center(single_peptite) - _packing_center(face_pep)
        face_pep = _translate_peptide(face_pep, face_delta[0], 0.0, face_delta[2])
        face_pep = _register_parallel_ca_yz(
            face_pep, single_peptite, residue_offset=1, direction=registry_direction
        )
    else:
        face_pep = None

    sheet1_parts = []
    for chain in range(chains_per_sheet):
        if anti_pep is not None and chain % 2 == 1:
            new_peptide = anti_pep.copy()
        elif face_pep is not None and chain % 2 == 1:
            new_peptide = face_pep.copy()
        else:
            new_peptide = single_peptite.copy()
        new_peptide["x"] = new_peptide["x"] + chain * strand_dist
        new_peptide["anum"] = new_peptide["anum"] + chain * anum_max
        new_peptide["resid"] = new_peptide["resid"] + chain * resid_max
        sheet1_parts.append(new_peptide)

    sheet1 = pd.concat(sheet1_parts, ignore_index=True)
    sheet2 = sheet1.copy()
    if class_axis[classes] != "none":
        sheet2 = rotation_coordinates(sheet2, class_axis[classes], np.pi)
    sheet2 = _align_strand_xy(
        sheet2, sheet1, same_direction=None, chain_length=chain_length
    )

    sheet1_center = _chain_center(sheet1, chain_length, chain_index=0)
    sheet2_center = _chain_center(sheet2, chain_length, chain_index=0)
    dx = (
        _packing_axis_midpoint(sheet1, "x")
        - _packing_axis_midpoint(sheet2, "x")
        + sheet_x_shift
    )
    dz = -sheet_dist if (core_residues == "o" and classes in [1, 3, 5]) else sheet_dist
    sheet2 = _translate_peptide(
        sheet2,
        dx,
        sheet1_center[1] + sheet_y_shift - sheet2_center[1],
        sheet1_center[2] + dz - sheet2_center[2],
    )
    sheet2["anum"] = sheet2["anum"] + chains_per_sheet * anum_max
    sheet2["resid"] = sheet2["resid"] + chains_per_sheet * resid_max

    peptides = pd.concat([sheet1, sheet2], ignore_index=True)
    peptides["anum"] = range(1, len(peptides) + 1)
    _write_packed_peptides(peptides, chain_length, format_flag, o_filename)
    return peptides


def _normalize_two_component_pattern(pattern, name):
    """Function:
        Clean and validate an A/B sheet pattern.

    Parameters
    ----------
    pattern : str
        Sequence of A and B letters that defines strand order in one sheet.
    name : str
        Parameter name used in pattern-validation error messages.

    Returns
    -------
    str
        String containing only upper-case A and B characters.
    """
    pattern = str(pattern).upper().replace(" ", "")
    if not pattern:
        raise ValueError(f"{name} cannot be empty")
    if any(letter not in ["A", "B"] for letter in pattern):
        raise ValueError(f"{name} can only contain A and B")
    return pattern


def _write_packed_peptides_two_component(peptides, chain_length, format_flag,
                                         o_filename):
    """Function:
        Write a packed two-component sheet and its A/B component labels to a PDB file.

    Parameters
    ----------
    peptides : pandas.DataFrame
        All packed peptide atom rows.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    format_flag : int
        0 selects PepAD format; 1 selects AMBER format.
    o_filename : str or pathlib.Path
        Output PDB file.

    Returns
    -------
    None
        Writes a two-component PDB file.
    """
    with open(o_filename, "w") as f:
        last_resid = None
        for _, row in peptides.iterrows():
            resid = int(row["resid"])
            aa_name = (
                row["aa_name"][1:]
                if format_flag == 1 and len(row["aa_name"]) == 4
                else row["aa_name"]
            )
            component = row["component"] if "component" in row else "A"
            if format_flag == 1 and last_resid is not None:
                if (last_resid % chain_length == 0) and (resid % chain_length == 1):
                    f.write("TER\n")
            output = "ATOM{:7d}{:^6}{:4}{:1}{:4d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                int(row["anum"]),
                row["atom_name"],
                aa_name,
                component,
                resid,
                float(row["x"]),
                float(row["y"]),
                float(row["z"]),
            )
            f.write(output)
            last_resid = resid


def _build_aligned_peptide_template(
    sequence: str, angles: Sequence[float],
    aligned_file: Union[str, Path], cap_flag: int,
    max_cb_tilt_deg: Optional[float] = None
) -> Tuple[int, pd.DataFrame]:
    """Function:
        Build, terminate, and align one peptide using an in-memory DataFrame.

    Parameters
    ----------
    sequence : str
        Peptide sequence using one-letter amino-acid codes.
    angles : list or tuple [phi, psi] in degrees
        Backbone angles for this template.
    aligned_file : str or pathlib.Path
        Persistent checking PDB for the aligned single strand.
    cap_flag : int
        0 uses uncapped termini and 1 or 2 selects a cap style.
    max_cb_tilt_deg : float or None
        Maximum allowed middle CA-CB tilt in degrees.

    Returns
    -------
    tuple
        Tuple (length, aligned_df), where length is int and aligned_df is
        pandas.DataFrame.
    """
    peptide_df = build_single_peptide_with_local_peptidebuilder(sequence, angles)
    if cap_flag == 0:
        peptide_df = add_NH_remove_OH(peptide_df)
        length = len(sequence)
    elif cap_flag in [1, 2]:
        peptide_df = add_caps(peptide_df, cap_flag)
        length = len(sequence) + 2
    else:
        raise ValueError("Cannot determine caps or not. Stop")

    aligned_df = alignment(
        length, 1, peptide_df, aligned_file, order=1,
        max_cb_tilt_deg=max_cb_tilt_deg
    )
    return length, aligned_df

def _two_component_variants(single_peptite, chain_length, classes,
                            registry_direction="minus"):
    """Function:
        Prepare the base, reversed, or face-flipped strand variants needed for one
        component.

    Parameters
    ----------
    single_peptite : pandas.DataFrame
        One aligned peptide strand used as a packing template.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    classes : int
        Steric-zipper class number from 1 to 10.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    dict
        Dictionary mapping variant names to pandas.DataFrame strand templates.
    """
    variants = {"base": single_peptite.copy()}
    if classes in [5, 6, 7, 8]:
        anti_axis = "z" if classes in [5, 6] else "x"
        anti_pep = rotation_coordinates(single_peptite, anti_axis, np.pi)
        anti_pep = _align_strand_xy(anti_pep, single_peptite, same_direction=False)
        if classes in [5, 6]:
            register_offset = 1 if chain_length % 2 == 0 else 0
        else:
            register_offset = 1 if chain_length % 2 == 1 else 0
        anti_pep = _register_antiparallel_ca_y(
            anti_pep,
            single_peptite,
            residue_offset=register_offset,
            direction=registry_direction,
        )
        center_delta = _packing_center(single_peptite) - _packing_center(anti_pep)
        variants["anti"] = _translate_peptide(
            anti_pep, center_delta[0], 0.0, center_delta[2]
        )
    if classes in [9, 10]:
        face_pep = rotation_coordinates(single_peptite, "y", np.pi)
        face_pep = _align_strand_xy(face_pep, single_peptite, same_direction=True)
        face_delta = _packing_center(single_peptite) - _packing_center(face_pep)
        face_pep = _translate_peptide(face_pep, face_delta[0], 0.0, face_delta[2])
        variants["face"] = _register_parallel_ca_yz(
            face_pep, single_peptite, residue_offset=1, direction=registry_direction
        )
    return variants


def _make_two_component_sheet(pattern, component_variants, classes, strand_dist,
                              anum_max, resid_max):
    """Function:
        Construct one sheet by placing peptide A or B according to an A/B pattern.

    Parameters
    ----------
    pattern : str
        Sequence of A and B letters that defines strand order in one sheet.
    component_variants : dict
        Maps component A/B and variant names to pandas.DataFrame strand templates.
    classes : int
        Steric-zipper class number from 1 to 10.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    anum_max : int
        Atom-number offset reserved for each added chain.
    resid_max : int
        Residue-number offset reserved for each added chain.

    Returns
    -------
    pandas.DataFrame
        Atom table containing one two-component sheet.
    """
    sheet_parts = []
    for chain, component in enumerate(pattern):
        variants = component_variants[component]
        if classes in [5, 6, 7, 8] and chain % 2 == 1:
            new_peptide = variants["anti"].copy()
        elif classes in [9, 10] and chain % 2 == 1:
            new_peptide = variants["face"].copy()
        else:
            new_peptide = variants["base"].copy()
        new_peptide["component"] = component
        new_peptide["x"] = new_peptide["x"] + chain * strand_dist
        new_peptide["anum"] = new_peptide["anum"] + chain * anum_max
        new_peptide["resid"] = new_peptide["resid"] + chain * resid_max
        sheet_parts.append(new_peptide)
    return pd.concat(sheet_parts, ignore_index=True)


def packpep_geometric_two_component(single_pep_A, single_pep_B, chain_length, classes,
                                    shifts, sheet1_pattern, sheet2_pattern, format_flag,
                                    o_filename, core_residues="e", strand_dist=4.8,
                                    sheet_dist=11.5, sheet_shift=0.0,
                                    registry_direction="minus"):
    """Function:
        Pack peptide A and peptide B into two uniform parallel or antiparallel sheets
        using two A/B patterns.

    Parameters
    ----------
    single_pep_A : pandas.DataFrame
        Aligned peptide-A strand template.
    single_pep_B : pandas.DataFrame
        Aligned peptide-B strand template.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    classes : int
        Steric-zipper class number from 1 to 10.
    shifts : int or float
        Move sheet 2 along y by this number of residue spacings.
    sheet1_pattern : str
        A/B strand pattern for sheet 1.
    sheet2_pattern : str
        A/B strand pattern for sheet 2.
    format_flag : int
        0 selects PepAD format; 1 selects AMBER format.
    o_filename : str or pathlib.Path
        Output PDB file.
    core_residues : str
        'e' packs even-numbered residues inside and 'o' packs odd-numbered residues
        inside.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    sheet_dist : float
        Sheet-to-sheet distance along z in Angstrom.
    sheet_shift : float
        Move sheet 2 along x in units of half the strand distance.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    pandas.DataFrame containing both packed two-component sheets
        Also writes o_filename.
    """
    class_axis = {
        1: "x",
        2: "none",
        3: "y",
        4: "z",
        5: "x",
        6: "none",
        7: "none",
        8: "z",
        9: "none",
        10: "x",
    }
    if classes not in class_axis:
        raise ValueError("classes must be an integer from 1 to 10")
    sheet1_pattern = _normalize_two_component_pattern(sheet1_pattern, "sheet1_pattern")
    sheet2_pattern = _normalize_two_component_pattern(sheet2_pattern, "sheet2_pattern")
    if len(sheet1_pattern) != len(sheet2_pattern):
        raise ValueError("sheet1_pattern and sheet2_pattern must have the same length")
    if len(_ca_trace(single_pep_A)) != len(_ca_trace(single_pep_B)):
        raise ValueError("Peptide A and peptide B must have the same residue length")
    if registry_direction not in ["minus", "plus"]:
        raise ValueError("registry_direction must be 'minus' or 'plus'")

    single_pep_B = _register_parallel_ca_yz(
        single_pep_B, single_pep_A, residue_offset=0
    )
    anum_max = max(int(single_pep_A["anum"].max()), int(single_pep_B["anum"].max()))
    resid_max = max(int(single_pep_A["resid"].max()), int(single_pep_B["resid"].max()))
    component_variants = {
        "A": _two_component_variants(
            single_pep_A, chain_length, classes, registry_direction
        ),
        "B": _two_component_variants(
            single_pep_B, chain_length, classes, registry_direction
        ),
    }

    sheet1 = _make_two_component_sheet(
        sheet1_pattern, component_variants, classes, strand_dist, anum_max, resid_max
    )
    sheet2 = _make_two_component_sheet(
        sheet2_pattern, component_variants, classes, strand_dist, anum_max, resid_max
    )
    if class_axis[classes] != "none":
        sheet2 = rotation_coordinates(sheet2, class_axis[classes], np.pi)
    sheet2 = _align_strand_xy(
        sheet2, sheet1, same_direction=None, chain_length=chain_length
    )

    residue_spacing = _estimate_residue_spacing(single_pep_A)
    sheet_y_shift = float(shifts) * residue_spacing
    sheet_x_shift = float(sheet_shift) * float(strand_dist) * 0.5
    sheet1_center = _chain_center(sheet1, chain_length, chain_index=0)
    sheet2_center = _chain_center(sheet2, chain_length, chain_index=0)
    dx = (
        _packing_axis_midpoint(sheet1, "x")
        - _packing_axis_midpoint(sheet2, "x")
        + sheet_x_shift
    )
    dz = -sheet_dist if (core_residues == "o" and classes in [1, 3, 5]) else sheet_dist
    sheet2 = _translate_peptide(
        sheet2,
        dx,
        sheet1_center[1] + sheet_y_shift - sheet2_center[1],
        sheet1_center[2] + dz - sheet2_center[2],
    )
    sheet2["anum"] = sheet2["anum"] + len(sheet1_pattern) * anum_max
    sheet2["resid"] = sheet2["resid"] + len(sheet1_pattern) * resid_max

    peptides = pd.concat([sheet1, sheet2], ignore_index=True)
    peptides["anum"] = range(1, len(peptides) + 1)
    _write_packed_peptides_two_component(
        peptides, chain_length, format_flag, o_filename
    )
    return peptides


def _make_two_component_hybrid_sheet(pattern, component_templates, strand_dist,
                                     anum_max, resid_max, alternate=False):
    """Function:
        Construct one hybrid sheet from A/B component templates and an A/B pattern.

    Parameters
    ----------
    pattern : str
        Sequence of A and B letters that defines strand order in one sheet.
    component_templates : dict
        Maps component A/B and template roles to pandas.DataFrame strands.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    anum_max : int
        Atom-number offset reserved for each added chain.
    resid_max : int
        Residue-number offset reserved for each added chain.
    alternate : bool
        When True, alternating chains use the secondary template.

    Returns
    -------
    pandas.DataFrame
        Atom table containing one hybrid component sheet.
    """
    sheet_parts = []
    for chain, component in enumerate(pattern):
        template_name = "secondary" if alternate and chain % 2 == 1 else "primary"
        new_peptide = component_templates[component][template_name].copy()
        new_peptide["component"] = component
        new_peptide["x"] = new_peptide["x"] + chain * strand_dist
        new_peptide["anum"] = new_peptide["anum"] + chain * anum_max
        new_peptide["resid"] = new_peptide["resid"] + chain * resid_max
        sheet_parts.append(new_peptide)
    return pd.concat(sheet_parts, ignore_index=True)


def packpep_parallel_antiparallel_two_component(single_para_A, single_para_B,
                                                single_anti_A, single_anti_B,
                                                chain_length, hybrid_type, shifts,
                                                sheet1_pattern, sheet2_pattern,
                                                format_flag, o_filename,
                                                core_residues="e", strand_dist=4.8,
                                                sheet_dist=11.5, sheet_shift=0.0,
                                                registry_direction="minus"):
    """Function:
        Pack two peptide components into one parallel sheet and one antiparallel sheet.

    Parameters
    ----------
    single_para_A : pandas.DataFrame
        Peptide-A template built with parallel angles.
    single_para_B : pandas.DataFrame
        Peptide-B template built with parallel angles.
    single_anti_A : pandas.DataFrame
        Peptide-A template built with antiparallel angles.
    single_anti_B : pandas.DataFrame
        Peptide-B template built with antiparallel angles.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    hybrid_type : int or None
        Hybrid packing type from 1 to 6, or None for a uniform sheet.
    shifts : int or float
        Move sheet 2 along y by this number of residue spacings.
    sheet1_pattern : str
        A/B strand pattern for sheet 1.
    sheet2_pattern : str
        A/B strand pattern for sheet 2.
    format_flag : int
        0 selects PepAD format; 1 selects AMBER format.
    o_filename : str or pathlib.Path
        Output PDB file.
    core_residues : str
        'e' packs even-numbered residues inside and 'o' packs odd-numbered residues
        inside.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    sheet_dist : float
        Sheet-to-sheet distance along z in Angstrom.
    sheet_shift : float
        Move sheet 2 along x in units of half the strand distance.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    pandas.DataFrame containing the packed two-component hybrid
        Also writes o_filename.
    """
    if hybrid_type not in [1, 2, 3, 4, 5, 6]:
        raise ValueError("hybrid_type must be an integer from 1 to 6")
    sheet1_pattern = _normalize_two_component_pattern(sheet1_pattern, "sheet1_pattern")
    sheet2_pattern = _normalize_two_component_pattern(sheet2_pattern, "sheet2_pattern")
    if len(sheet1_pattern) != len(sheet2_pattern):
        raise ValueError("sheet1_pattern and sheet2_pattern must have the same length")
    if registry_direction not in ["minus", "plus"]:
        raise ValueError("registry_direction must be 'minus' or 'plus'")

    templates = [single_para_A, single_para_B, single_anti_A, single_anti_B]
    ca_counts = [len(_ca_trace(template)) for template in templates]
    if len(set(ca_counts)) != 1:
        raise ValueError(
            "All two-component hybrid templates must have the same residue length"
        )

    single_para_B = _register_parallel_ca_yz(
        single_para_B, single_para_A, residue_offset=0
    )
    single_anti_B = _register_parallel_ca_yz(
        single_anti_B, single_anti_A, residue_offset=0
    )
    anum_max = max(int(template["anum"].max()) for template in templates)
    resid_max = max(int(template["resid"].max()) for template in templates)
    residue_spacing = _estimate_residue_spacing(single_para_A)
    sheet_y_shift = float(shifts) * residue_spacing
    sheet_x_shift = float(sheet_shift) * float(strand_dist) * 0.5

    sheet1_templates = {
        1: (False, False),
        2: (False, False),
        3: (False, False),
        4: (False, True),
        5: (True, False),
        6: (True, True),
    }
    sheet1_face_flipped, sheet1_direction_reversed = sheet1_templates[hybrid_type]
    sheet1_component_templates = {}
    for component, template in [("A", single_para_A), ("B", single_para_B)]:
        primary = _make_hybrid_strand_template(
            template,
            face_flipped=sheet1_face_flipped,
            direction_reversed=sheet1_direction_reversed,
        )
        sheet1_component_templates[component] = {
            "primary": primary,
            "secondary": primary.copy(),
        }

    sheet2_component_templates = {}
    for component, template in [("A", single_anti_A), ("B", single_anti_B)]:
        if hybrid_type == 1:
            primary = _make_hybrid_strand_template(
                template, face_flipped=True, direction_reversed=False
            )
        elif hybrid_type == 2:
            primary = template.copy()
        else:
            primary = _make_hybrid_strand_template(
                template, face_flipped=True, direction_reversed=False
            )
        if hybrid_type in [1, 2]:
            secondary = _class56_antiparallel_neighbor(
                primary, chain_length, registry_direction=registry_direction
            )
        else:
            secondary = _class78_antiparallel_neighbor(
                primary, chain_length, registry_direction=registry_direction
            )
        sheet2_component_templates[component] = {
            "primary": primary,
            "secondary": secondary,
        }

    sheet1 = _make_two_component_hybrid_sheet(
        sheet1_pattern,
        sheet1_component_templates,
        strand_dist,
        anum_max,
        resid_max,
        alternate=False,
    )
    sheet2 = _make_two_component_hybrid_sheet(
        sheet2_pattern,
        sheet2_component_templates,
        strand_dist,
        anum_max,
        resid_max,
        alternate=True,
    )

    sheet1_center = _chain_center(sheet1, chain_length, chain_index=0)
    sheet2_center = _chain_center(sheet2, chain_length, chain_index=0)
    dx = (
        _packing_axis_midpoint(sheet1, "x")
        - _packing_axis_midpoint(sheet2, "x")
        + sheet_x_shift
    )
    dz = -sheet_dist if core_residues == "o" else sheet_dist
    sheet2 = _translate_peptide(
        sheet2,
        dx,
        sheet1_center[1] + sheet_y_shift - sheet2_center[1],
        sheet1_center[2] + dz - sheet2_center[2],
    )
    sheet2["anum"] = sheet2["anum"] + len(sheet1_pattern) * anum_max
    sheet2["resid"] = sheet2["resid"] + len(sheet1_pattern) * resid_max

    peptides = pd.concat([sheet1, sheet2], ignore_index=True)
    peptides["anum"] = range(1, len(peptides) + 1)
    _write_packed_peptides_two_component(
        peptides, chain_length, format_flag, o_filename
    )
    return peptides


def _register_parallel_ca_y(moving_df, target_df):
    """Function:
        Translate one parallel strand along y so its CA positions match a target strand.

    Parameters
    ----------
    moving_df : pandas.DataFrame
        Strand or sheet coordinates that will be moved.
    target_df : pandas.DataFrame
        Fixed reference strand or sheet coordinates.

    Returns
    -------
    pandas.DataFrame
        Peptide atom table translated to match the target CA y positions.
    """
    moving_y = _ca_trace(moving_df)["y"].to_numpy(dtype=np.float64)
    target_y = _ca_trace(target_df)["y"].to_numpy(dtype=np.float64)
    n = min(len(moving_y), len(target_y))
    if n == 0:
        return moving_df
    return _translate_peptide(moving_df, dy=float(np.mean(target_y[:n] - moving_y[:n])))


def _class56_antiparallel_neighbor(single_peptite, chain_length,
                                   registry_direction="minus"):
    """Function:
        Create the neighboring antiparallel strand and registry required for Class 5 or 6.

    Parameters
    ----------
    single_peptite : pandas.DataFrame
        One aligned peptide strand used as a packing template.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the Class 5/6 neighboring strand.
    """
    anti_pep = rotation_coordinates(single_peptite, "z", np.pi)
    anti_pep = _align_strand_xy(anti_pep, single_peptite, same_direction=False)
    register_offset = 1 if chain_length % 2 == 0 else 0
    anti_pep = _register_antiparallel_ca_y(
        anti_pep,
        single_peptite,
        residue_offset=register_offset,
        direction=registry_direction,
    )
    center_delta = _packing_center(single_peptite) - _packing_center(anti_pep)
    return _translate_peptide(anti_pep, center_delta[0], 0.0, center_delta[2])


def _class78_antiparallel_neighbor(single_peptite, chain_length,
                                   registry_direction="minus"):
    """Function:
        Create the neighboring antiparallel strand and registry required for Class 7 or
        8.

    Parameters
    ----------
    single_peptite : pandas.DataFrame
        One aligned peptide strand used as a packing template.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the Class 7/8 neighboring strand.
    """
    anti_pep = rotation_coordinates(single_peptite, "x", np.pi)
    anti_pep = _align_strand_xy(anti_pep, single_peptite, same_direction=False)
    register_offset = 1 if chain_length % 2 == 1 else 0
    anti_pep = _register_antiparallel_ca_y(
        anti_pep,
        single_peptite,
        residue_offset=register_offset,
        direction=registry_direction,
    )
    center_delta = _packing_center(single_peptite) - _packing_center(anti_pep)
    return _translate_peptide(anti_pep, center_delta[0], 0.0, center_delta[2])


def _make_hybrid_strand_template(single_peptite, face_flipped=False,
                                 direction_reversed=False):
    """Function:
        Flip the face or reverse the direction of one strand template for hybrid
        packing.

    Parameters
    ----------
    single_peptite : pandas.DataFrame
        One aligned peptide strand used as a packing template.
    face_flipped : bool
        When True, turn the side-chain face to the opposite side.
    direction_reversed : bool
        When True, reverse the strand N-to-C direction.

    Returns
    -------
    pandas.DataFrame
        Atom table containing the transformed strand template.
    """
    if face_flipped and direction_reversed:
        template = rotation_coordinates(single_peptite, "x", np.pi)
    elif face_flipped:
        template = rotation_coordinates(single_peptite, "y", np.pi)
    elif direction_reversed:
        template = rotation_coordinates(single_peptite, "z", np.pi)
    else:
        return single_peptite.copy()

    template = _align_strand_xy(
        template, single_peptite, same_direction=not direction_reversed
    )
    if direction_reversed:
        template = _register_antiparallel_ca_y(template, single_peptite)
    else:
        template = _register_parallel_ca_y(template, single_peptite)
    center_delta = _packing_center(single_peptite) - _packing_center(template)
    return _translate_peptide(template, center_delta[0], 0.0, center_delta[2])


def _make_sheet_from_templates(primary_pep, secondary_pep, chains_per_sheet,
                               strand_dist, anum_max, resid_max):
    """Function:
        Build one sheet by alternating primary and secondary strand templates along x.

    Parameters
    ----------
    primary_pep : pandas.DataFrame
        Primary strand template.
    secondary_pep : pandas.DataFrame
        Secondary strand template used on alternating chains.
    chains_per_sheet : int
        Number of peptide strands placed in each sheet.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    anum_max : int
        Atom-number offset reserved for each added chain.
    resid_max : int
        Residue-number offset reserved for each added chain.

    Returns
    -------
    pandas.DataFrame
        Atom table containing one packed sheet.
    """
    sheet_parts = []
    for chain in range(chains_per_sheet):
        new_peptide = primary_pep.copy() if chain % 2 == 0 else secondary_pep.copy()
        new_peptide["x"] = new_peptide["x"] + chain * strand_dist
        new_peptide["anum"] = new_peptide["anum"] + chain * anum_max
        new_peptide["resid"] = new_peptide["resid"] + chain * resid_max
        sheet_parts.append(new_peptide)
    return pd.concat(sheet_parts, ignore_index=True)


def packpep_parallel_antiparallel(single_pep_para, single_pep_anti=None,
                                  chain_length=None, hybrid_type=None, shifts=None,
                                  chains_per_sheet=None, format_flag=None,
                                  o_filename=None, core_residues="e", strand_dist=4.8,
                                  sheet_dist=11.5, sheet_shift=0.0,
                                  registry_direction="minus"):
    """Function:
        Pack one parallel-angle template and one antiparallel-angle template into a
        hybrid two-sheet structure.

    Parameters
    ----------
    single_pep_para : pandas.DataFrame
        Peptide template built with parallel beta-strand angles.
    single_pep_anti : pandas.DataFrame
        Peptide template built with antiparallel beta-strand angles.
    chain_length : int
        Number of residue records assigned to one packed chain, including caps when
        present.
    hybrid_type : int or None
        Hybrid packing type from 1 to 6, or None for a uniform sheet.
    shifts : int or float
        Move sheet 2 along y by this number of residue spacings.
    chains_per_sheet : int
        Number of peptide strands placed in each sheet.
    format_flag : int
        0 selects PepAD format; 1 selects AMBER format.
    o_filename : str or pathlib.Path
        Output PDB file.
    core_residues : str
        'e' packs even-numbered residues inside and 'o' packs odd-numbered residues
        inside.
    strand_dist : float
        Neighboring-strand distance along x in Angstrom.
    sheet_dist : float
        Sheet-to-sheet distance along z in Angstrom.
    sheet_shift : float
        Move sheet 2 along x in units of half the strand distance.
    registry_direction : str
        'minus' or 'plus' selects the y direction of a one-residue registry shift.

    Returns
    -------
    pandas.DataFrame containing the packed hybrid sheets
        Also writes o_filename.
    """
    if not isinstance(single_pep_anti, pd.DataFrame):
        old_chain_length = single_pep_anti
        old_hybrid_type = chain_length
        old_shifts = hybrid_type
        old_chains_per_sheet = shifts
        old_format_flag = chains_per_sheet
        old_o_filename = format_flag
        old_core_residues = o_filename
        single_pep_anti = single_pep_para
        chain_length = old_chain_length
        hybrid_type = old_hybrid_type
        shifts = old_shifts
        chains_per_sheet = old_chains_per_sheet
        format_flag = old_format_flag
        o_filename = old_o_filename
        if old_core_residues is not None:
            core_residues = old_core_residues

    if hybrid_type not in [1, 2, 3, 4, 5, 6]:
        raise ValueError("hybrid_type must be an integer from 1 to 6")
    if registry_direction not in ["minus", "plus"]:
        raise ValueError("registry_direction must be 'minus' or 'plus'")
    if len(_ca_trace(single_pep_para)) != len(_ca_trace(single_pep_anti)):
        raise ValueError(
            "Parallel and antiparallel templates must have the same number of residues"
        )

    anum_max = max(
        int(single_pep_para["anum"].max()), int(single_pep_anti["anum"].max())
    )
    resid_max = max(
        int(single_pep_para["resid"].max()), int(single_pep_anti["resid"].max())
    )
    residue_spacing = _estimate_residue_spacing(single_pep_para)
    sheet_y_shift = float(shifts) * residue_spacing
    sheet_x_shift = float(sheet_shift) * float(strand_dist) * 0.5

    sheet1_templates = {
        1: (False, False),
        2: (False, False),
        3: (False, False),
        4: (False, True),
        5: (True, False),
        6: (True, True),
    }
    sheet1_face_flipped, sheet1_direction_reversed = sheet1_templates[hybrid_type]
    sheet1_primary = _make_hybrid_strand_template(
        single_pep_para,
        face_flipped=sheet1_face_flipped,
        direction_reversed=sheet1_direction_reversed,
    )
    sheet1_secondary = sheet1_primary.copy()
    sheet1 = _make_sheet_from_templates(
        sheet1_primary,
        sheet1_secondary,
        chains_per_sheet,
        strand_dist,
        anum_max,
        resid_max,
    )

    if hybrid_type == 1:
        sheet2_primary = _make_hybrid_strand_template(
            single_pep_anti, face_flipped=True, direction_reversed=False
        )
    else:
        sheet2_primary = (
            single_pep_anti.copy()
            if hybrid_type == 2
            else _make_hybrid_strand_template(
                single_pep_anti, face_flipped=True, direction_reversed=False
            )
        )

    if hybrid_type in [1, 2]:
        sheet2_secondary = _class56_antiparallel_neighbor(
            sheet2_primary, chain_length, registry_direction=registry_direction
        )
    else:
        sheet2_secondary = _class78_antiparallel_neighbor(
            sheet2_primary, chain_length, registry_direction=registry_direction
        )
    sheet2 = _make_sheet_from_templates(
        sheet2_primary,
        sheet2_secondary,
        chains_per_sheet,
        strand_dist,
        anum_max,
        resid_max,
    )

    sheet1_center = _chain_center(sheet1, chain_length, chain_index=0)
    sheet2_center = _chain_center(sheet2, chain_length, chain_index=0)
    dx = (
        _packing_axis_midpoint(sheet1, "x")
        - _packing_axis_midpoint(sheet2, "x")
        + sheet_x_shift
    )
    dz = -sheet_dist if core_residues == "o" else sheet_dist
    sheet2 = _translate_peptide(
        sheet2,
        dx,
        sheet1_center[1] + sheet_y_shift - sheet2_center[1],
        sheet1_center[2] + dz - sheet2_center[2],
    )
    sheet2["anum"] = sheet2["anum"] + chains_per_sheet * anum_max
    sheet2["resid"] = sheet2["resid"] + chains_per_sheet * resid_max

    peptides = pd.concat([sheet1, sheet2], ignore_index=True)
    peptides["anum"] = range(1, len(peptides) + 1)
    _write_packed_peptides(peptides, chain_length, format_flag, o_filename)
    return peptides


# %% Cell v122-input-help-system
SHEET_INPUT_SCHEMA = {
    "sequence": {
        "type": "str",
        "default": "required",
        "description": "Only peptide sequence, or peptide A in a two-sequence build.",
    },
    "output_name": {
        "type": "str",
        "default": "peptide_sheets",
        "description": "Output PDB filename. The .pdb suffix is added when absent.",
    },
    "sequence_b": {
        "type": "str or None",
        "default": None,
        "description": "Optional peptide B. It must have the same length as peptide A.",
    },
    "pattern1": {
        "type": "str or None",
        "default": None,
        "description": "Sheet-1 A/B strand pattern. Required with sequence_b.",
    },
    "pattern2": {
        "type": "str or None",
        "default": None,
        "description": "Sheet-2 A/B strand pattern. Required with sequence_b.",
    },
    "classes": {
        "type": "int or None",
        "default": None,
        "choices": tuple(range(1, 11)),
        "description": "Steric-zipper class. Do not combine with hybrid_type.",
    },
    "hybrid_type": {
        "type": "int or None",
        "default": None,
        "choices": tuple(range(1, 7)),
        "description": "Parallel-antiparallel type. Do not combine with classes.",
    },
    "core_residues": {
        "type": "{'e', 'o'}",
        "default": "e",
        "choices": ("e", "o"),
        "description": "Residue parity packed inside: 'e' for even or 'o' for odd.",
    },
    "chains": {
        "type": "int",
        "default": 8,
        "description": "Strands per sheet for one sequence. Patterns replace this.",
    },
    "format_flag": {
        "type": "{0, 1}",
        "default": 0,
        "choices": (0, 1),
        "description": "0 writes PepAD format; 1 writes AMBER format.",
    },
    "cap_flag": {
        "type": "{0, 1, 2}",
        "default": 0,
        "choices": (0, 1, 2),
        "description": "0 is uncapped; 1 adds ACE+NME; 2 adds ACE+NHE.",
    },
    "y_shift": {
        "type": "float",
        "default": 0.0,
        "description": "Sheet-2 registry shift along y, in residue-spacing units.",
    },
    "registry_direction": {
        "type": "{'minus', 'plus'}",
        "default": "minus",
        "choices": ("minus", "plus"),
        "description": "Direction used when a one-residue registry shift is required.",
    },
    "x_spacing": {
        "type": "float",
        "default": 4.8,
        "description": "Neighboring-strand spacing along x, in Angstrom.",
    },
    "z_spacing": {
        "type": "float",
        "default": 11.5,
        "description": "Sheet-to-sheet spacing along z, in Angstrom.",
    },
    "x_shift": {
        "type": "float",
        "default": 0.0,
        "description": "Extra sheet-2 offset along x, in half-x_spacing units.",
    },
    "max_cb_tilt_deg": {
        "type": "float or None",
        "default": 10.0,
        "description": (
            "Middle CA-CB tilt limit for parallel-angle strands. "
            "None disables the limit."
        ),
    },
}


UNIFORM_CLASS_HELP = {
    1: "parallel strands; sheet 2 is rotated 180 degrees around x",
    2: "parallel strands; sheet 2 keeps the same xyz orientation",
    3: "parallel strands; sheet 2 is rotated 180 degrees around y",
    4: "parallel strands; sheet 2 is rotated 180 degrees around z",
    5: (
        "antiparallel neighboring strands with Class 5/6 registry; "
        "sheet 2 rotates around x"
    ),
    6: (
        "antiparallel neighboring strands with Class 5/6 registry; "
        "sheet 2 keeps its orientation"
    ),
    7: (
        "antiparallel neighboring strands with Class 7/8 registry; "
        "sheet 2 keeps its orientation"
    ),
    8: (
        "antiparallel neighboring strands with Class 7/8 registry; "
        "sheet 2 rotates around z"
    ),
    9: (
        "parallel neighboring strands with alternating faces; "
        "sheet 2 keeps its orientation"
    ),
    10: "parallel neighboring strands with alternating faces; sheet 2 rotates around x",
}


HYBRID_TYPE_HELP = {
    1: "standard parallel sheet + face-flipped Class 5/6-style antiparallel sheet",
    2: "standard parallel sheet + unflipped Class 5/6-style antiparallel sheet",
    3: "standard parallel sheet + face-flipped Class 7/8-style antiparallel sheet",
    4: (
        "direction-reversed parallel sheet + face-flipped "
        "Class 7/8-style antiparallel sheet"
    ),
    5: "face-flipped parallel sheet + face-flipped Class 7/8-style antiparallel sheet",
    6: (
        "face-flipped and direction-reversed parallel sheet + "
        "face-flipped Class 7/8-style antiparallel sheet"
    ),
}


def get_sheet_builder_help():
    """Function:
        Create complete help text for one- or two-sequence sheet building.

    Returns
    -------
    str
        Parameter descriptions, packing-type descriptions, and usage examples.
    """
    lines = [
        "Initial Structure Builder v1.23",
        "",
        "Supported builds:",
        "  1. One sequence with a steric-zipper class from 1 to 10.",
        "  2. One sequence with a parallel-antiparallel type from 1 to 6.",
        "  3. Two sequences with A/B patterns and a class from 1 to 10.",
        "  4. Two sequences with A/B patterns and a hybrid type from 1 to 6.",
        "",
        "Geometry selection:",
        "  Supply either classes or hybrid_type, but not both.",
        "  A direct notebook call uses Class 1 when both are omitted.",
        "  The terminal requires -c/--class or -H/--hybrid.",
        "",
        "Parameters:",
    ]
    for name, info in SHEET_INPUT_SCHEMA.items():
        choices = info.get("choices")
        choice_text = "" if choices is None else f"; choices={list(choices)}"
        lines.append(
            f"  {name} [{info['type']}; default={info['default']}{choice_text}]"
        )
        lines.append(f"      {info['description']}")

    lines.extend(["", "Steric-zipper classes:"])
    lines.extend(
        f"  {number}: {description}"
        for number, description in UNIFORM_CLASS_HELP.items()
    )
    lines.extend(["", "Parallel-antiparallel hybrid types:"])
    lines.extend(
        f"  {number}: {description}"
        for number, description in HYBRID_TYPE_HELP.items()
    )

    lines.extend(
        [
            "",
            "Notebook help:",
            "  show_sheet_builder_help()",
            "  show_input_help()",
            "  show_input_help('y_shift')",
            "",
            "Notebook builds:",
            "  build_sheets('GNNQQNY', 'class1', classes=1, chains=4)",
            "  build_sheets('GNNQQNY', 'hybrid1', hybrid_type=1, chains=4)",
            "  build_sheets('GNNQQNY', 'two_class', sequence_b='AAAAAAA',",
            "               pattern1='ABAB', pattern2='ABAB', classes=1)",
            "  build_sheets('GNNQQNY', 'two_hybrid', sequence_b='AAAAAAA',",
            "               pattern1='AAAA', pattern2='BBBB', hybrid_type=1)",
            "",
            "Terminal help after exporting to Python:",
            "  python initial_structure_builder.py -h",
            "  python initial_structure_builder.py -s GNNQQNY -c 1",
            "  python initial_structure_builder.py -s GNNQQNY -H 1",
            "  python initial_structure_builder.py -s GNNQQNY -b AAAAAAA",
            "      -p1 ABAB -p2 ABAB -c 1",
            "  python initial_structure_builder.py -s GNNQQNY -b AAAAAAA",
            "      -p1 AAAA -p2 BBBB -H 1",
        ]
    )
    return "\n".join(lines)


def show_sheet_builder_help():
    """Function:
        Print complete help for one- or two-sequence sheet building.

    Returns
    -------
    None
        The help text is printed in the notebook.
    """
    print(get_sheet_builder_help())


def show_input_help(parameter=None):
    """Function:
        Print help for one sheet-builder input.

    Parameters
    ----------
    parameter : str or None, default=None
        Input name, such as 'y_shift' or 'x_spacing'. None lists all input names.

    Returns
    -------
    None
        The requested input help is printed in the notebook.
    """
    if parameter is None:
        print("Available inputs:")
        print("  " + "\n  ".join(SHEET_INPUT_SCHEMA))
        print("\nUse show_input_help('input_name') for one detailed description.")
        return

    parameter = str(parameter).strip()
    if parameter not in SHEET_INPUT_SCHEMA:
        available = ", ".join(SHEET_INPUT_SCHEMA)
        raise ValueError(f"Unknown input '{parameter}'. Available inputs: {available}")

    info = SHEET_INPUT_SCHEMA[parameter]
    print(parameter)
    print(f"  type: {info['type']}")
    print(f"  default: {info['default']}")
    if "choices" in info:
        print(f"  choices: {list(info['choices'])}")
    print(f"  meaning: {info['description']}")


def _finite_float(value: Union[int, float], name: str) -> float:
    """Function:
        Convert one numeric builder input to a finite float.

    Parameters
    ----------
    value : int or float
        Numeric value to convert.
    name : str
        Parameter name used in an error message.

    Returns
    -------
    float
        Validated finite value.
    """
    try:
        converted = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a number") from exc
    if not np.isfinite(converted):
        raise ValueError(f"{name} must be finite")
    return converted


def validate_sheet_inputs(
    sequence: str, output_name: Union[str, Path] = "peptide_sheets", *,
    mode: str = "uniform", classes: Optional[int] = 1,
    hybrid_type: Optional[int] = None, core_residues: str = "e",
    chains: int = 8, format_flag: int = 0, cap_flag: int = 0,
    y_shift: float = 0.0, registry_direction: str = "minus",
    x_spacing: float = 4.8, z_spacing: float = 11.5,
    x_shift: float = 0.0, max_cb_tilt_deg: Optional[float] = 10.0
) -> dict:
    """Function:
        Validate and normalize inputs for sheet building.

    Parameters
    ----------
    sequence : str
        Peptide sequence in one-letter amino-acid codes.
    output_name : str, default='peptide_sheets'
        Output PDB filename.
    mode : {'uniform', 'hybrid'}, default='uniform'
        Builder execution mode.
    classes : int, default=1
        Uniform class from 1 to 10. Ignored when mode='hybrid'.
    hybrid_type : int or None, default=None
        Hybrid classification from 1 to 6.
    core_residues : {'e', 'o'}, default='e'
        Residue parity packed inside: 'e' for even or 'o' for odd.
    chains : int, default=8
        Number of strands per sheet.
    format_flag : {0, 1}, default=0
        0 writes PepAD format; 1 writes AMBER format.
    cap_flag : {0, 1, 2}, default=0
        0 is uncapped; 1 adds ACE+NME; 2 adds ACE+NHE.
    y_shift : float, default=0.0
        Sheet-2 registry shift along y, in residue-spacing units.
    registry_direction : {'minus', 'plus'}, default='minus'
        Direction used when a one-residue registry shift is required.
    x_spacing : float, default=4.8
        Neighboring-strand spacing along x, in Angstrom.
    z_spacing : float, default=11.5
        Sheet-to-sheet spacing along z, in Angstrom.
    x_shift : float, default=0.0
        Extra sheet-2 offset along x, in half-x_spacing units.
    max_cb_tilt_deg : float or None, default=10.0
        Middle CA-CB tilt limit for parallel-angle strands.

    Returns
    -------
    dict
        Normalized values ready for build_sheets().
    """
    mode = str(mode).lower().strip()
    if mode not in ("uniform", "hybrid"):
        raise ValueError("mode must be 'uniform' or 'hybrid'")

    sequence = "".join(str(sequence).split()).upper()
    if not sequence:
        raise ValueError("sequence cannot be empty")
    allowed = set("ARNDCQEGHILKMFPSTWYV")
    invalid = sorted(set(sequence) - allowed)
    if invalid:
        raise ValueError(
            f"sequence contains invalid one-letter amino-acid codes: {invalid}"
        )

    output_name = str(output_name).strip()
    if not output_name:
        raise ValueError("output_name cannot be empty")

    if (
        isinstance(chains, bool)
        or not isinstance(chains, (int, np.integer))
        or int(chains) < 1
    ):
        raise ValueError("chains must be a positive integer")
    chains = int(chains)

    if mode == "uniform":
        if isinstance(classes, bool) or not isinstance(classes, (int, np.integer)):
            raise ValueError("classes must be an integer from 1 to 10 in uniform mode")
        classes = int(classes)
        if classes not in UNIFORM_CLASS_HELP:
            raise ValueError("classes must be an integer from 1 to 10 in uniform mode")
        if hybrid_type is not None:
            raise ValueError("hybrid_type must be None in uniform mode")
    else:
        if isinstance(hybrid_type, bool) or not isinstance(
            hybrid_type, (int, np.integer)
        ):
            raise ValueError(
                "hybrid_type must be an integer from 1 to 6 in hybrid mode"
            )
        hybrid_type = int(hybrid_type)
        if hybrid_type not in HYBRID_TYPE_HELP:
            raise ValueError(
                "hybrid_type must be an integer from 1 to 6 in hybrid mode"
            )
        classes = None

    core_residues = str(core_residues).lower().strip()
    if core_residues not in ("e", "o"):
        raise ValueError("core_residues must be 'e' or 'o'")
    if format_flag not in (0, 1):
        raise ValueError("format_flag must be 0 or 1")
    if cap_flag not in (0, 1, 2):
        raise ValueError("cap_flag must be 0, 1, or 2")

    registry_direction = str(registry_direction).lower().strip()
    if registry_direction not in ("minus", "plus"):
        raise ValueError("registry_direction must be 'minus' or 'plus'")

    y_shift = _finite_float(y_shift, "y_shift")
    x_spacing = _finite_float(x_spacing, "x_spacing")
    z_spacing = _finite_float(z_spacing, "z_spacing")
    x_shift = _finite_float(x_shift, "x_shift")
    if x_spacing <= 0.0:
        raise ValueError("x_spacing must be larger than 0 Angstrom")
    if z_spacing <= 0.0:
        raise ValueError("z_spacing must be larger than 0 Angstrom")

    if max_cb_tilt_deg is not None:
        max_cb_tilt_deg = _finite_float(max_cb_tilt_deg, "max_cb_tilt_deg")
        if not 0.0 <= max_cb_tilt_deg <= 90.0:
            raise ValueError(
                "max_cb_tilt_deg must be between 0 and 90 degrees, or None"
            )

    return {
        "sequence": sequence,
        "output_name": output_name,
        "mode": mode,
        "classes": classes,
        "hybrid_type": hybrid_type,
        "core_residues": core_residues,
        "chains": chains,
        "format_flag": int(format_flag),
        "cap_flag": int(cap_flag),
        "y_shift": y_shift,
        "registry_direction": registry_direction,
        "x_spacing": x_spacing,
        "z_spacing": z_spacing,
        "x_shift": x_shift,
        "max_cb_tilt_deg": max_cb_tilt_deg,
    }


def build_sheets(
    sequence: str, output_name: Union[str, Path] = "peptide_sheets", *,
    sequence_b: Optional[str] = None, pattern1: Optional[str] = None,
    pattern2: Optional[str] = None, classes: Optional[int] = None,
    hybrid_type: Optional[int] = None, core_residues: str = "e",
    chains: int = 8, format_flag: int = 0, cap_flag: int = 0,
    y_shift: float = 0.0, registry_direction: str = "minus",
    x_spacing: float = 4.8, z_spacing: float = 11.5,
    x_shift: float = 0.0, max_cb_tilt_deg: Optional[float] = 10.0
) -> pd.DataFrame:
    """Function:
        Build one- or two-sequence uniform or hybrid peptide sheets.

    Parameters
    ----------
    sequence : str
        Only peptide sequence, or peptide A when sequence_b is supplied.
    output_name : str, default='peptide_sheets'
        Output PDB filename.
    sequence_b : str or None, default=None
        Optional peptide B sequence with the same residue length as peptide A.
    pattern1 : str or None, default=None
        Sheet-1 A/B pattern. Required when sequence_b is supplied.
    pattern2 : str or None, default=None
        Sheet-2 A/B pattern. Required when sequence_b is supplied.
    classes : int or None, default=None
        Steric-zipper class from 1 to 10. Class 1 is used when neither geometry
        selector is supplied.
    hybrid_type : int or None, default=None
        Parallel-antiparallel hybrid type from 1 to 6.
    core_residues : {'e', 'o'}, default='e'
        Residue parity packed inside.
    chains : int, default=8
        Number of strands per sheet for a one-sequence build. Pattern lengths
        set the strand counts for a two-sequence build.
    format_flag : {0, 1}, default=0
        PepAD or AMBER format.
    cap_flag : {0, 1, 2}, default=0
        Uncapped, ACE+NME, or ACE+NHE termini.
    y_shift : float, default=0.0
        Sheet-2 registry shift along y, in residue-spacing units.
    registry_direction : {'minus', 'plus'}, default='minus'
        Direction used for a required one-residue registry shift.
    x_spacing : float, default=4.8
        Neighboring-strand spacing along x, in Angstrom.
    z_spacing : float, default=11.5
        Sheet-to-sheet spacing along z, in Angstrom.
    x_shift : float, default=0.0
        Extra sheet-2 offset along x, in half-x_spacing units.
    max_cb_tilt_deg : float or None, default=10.0
        Middle CA-CB tilt limit for parallel-angle strands.

    Returns
    -------
    pandas.DataFrame
        Packed atoms from both sheets. The function also writes the output PDB file.
    """
    if classes is None and hybrid_type is None:
        classes = 1
    elif classes is not None and hybrid_type is not None:
        raise ValueError("Use either classes or hybrid_type, not both")

    if sequence_b and not (pattern1 and pattern2):
        raise ValueError("sequence_b requires both pattern1 and pattern2")
    if not sequence_b and (pattern1 or pattern2):
        raise ValueError("pattern1 and pattern2 require sequence_b")

    mode = "hybrid" if hybrid_type is not None else "uniform"
    validation_chains = 1 if sequence_b else chains
    config = validate_sheet_inputs(
        sequence,
        output_name, # after the second term, must enter inputs by names
        mode=mode,
        classes=classes,
        hybrid_type=hybrid_type,
        core_residues=core_residues,
        chains=validation_chains,
        format_flag=format_flag,
        cap_flag=cap_flag,
        y_shift=y_shift,
        registry_direction=registry_direction,
        x_spacing=x_spacing,
        z_spacing=z_spacing,
        x_shift=x_shift,
        max_cb_tilt_deg=max_cb_tilt_deg,
    )

    sequence_a = config["sequence"]
    if sequence_b:
        config_b = validate_sheet_inputs(
            sequence_b,
            output_name,
            mode=mode,
            classes=classes,
            hybrid_type=hybrid_type,
            core_residues=core_residues,
            chains=1,
            format_flag=format_flag,
            cap_flag=cap_flag,
            y_shift=y_shift,
            registry_direction=registry_direction,
            x_spacing=x_spacing,
            z_spacing=z_spacing,
            x_shift=x_shift,
            max_cb_tilt_deg=max_cb_tilt_deg,
        )
        sequence_b = config_b["sequence"]
        if len(sequence_a) != len(sequence_b):
            raise ValueError("Peptide A and peptide B must have equal lengths")
        pattern1 = _normalize_two_component_pattern(pattern1, "pattern1")
        pattern2 = _normalize_two_component_pattern(pattern2, "pattern2")

    output_file = config["output_name"]
    if not output_file.lower().endswith(".pdb"):
        output_file += ".pdb"

    if mode == "uniform":
        parallel_classes = (1, 2, 3, 4, 9, 10)
        angles = [-119, 113] if config["classes"] in parallel_classes else [-139, 135]
        tilt_limit = (
            config["max_cb_tilt_deg"] if config["classes"] in parallel_classes else None
        )
        template_name = "para" if config["classes"] in parallel_classes else "anti"
        length, peptide_a = _build_aligned_peptide_template(
            sequence_a, angles, f"{template_name}_A.pdb", config["cap_flag"],
            max_cb_tilt_deg=tilt_limit
        )

        if not sequence_b:
            result = packpep_geometric(
                peptide_a,
                length,
                config["classes"],
                config["y_shift"],
                config["chains"],
                config["format_flag"],
                output_file,
                core_residues=config["core_residues"],
                strand_dist=config["x_spacing"],
                sheet_dist=config["z_spacing"],
                sheet_shift=config["x_shift"],
                registry_direction=config["registry_direction"],
            )
        else:
            length_b, peptide_b = _build_aligned_peptide_template(
                sequence_b, angles, f"{template_name}_B.pdb",
                config["cap_flag"], max_cb_tilt_deg=tilt_limit
            )
            if length != length_b:
                raise ValueError("Peptide A and peptide B have different lengths")
            result = packpep_geometric_two_component(
                peptide_a,
                peptide_b,
                length,
                config["classes"],
                config["y_shift"],
                pattern1,
                pattern2,
                config["format_flag"],
                output_file,
                core_residues=config["core_residues"],
                strand_dist=config["x_spacing"],
                sheet_dist=config["z_spacing"],
                sheet_shift=config["x_shift"],
                registry_direction=config["registry_direction"],
            )
    else:
        length, parallel_a = _build_aligned_peptide_template(
            sequence_a, [-119, 113], "para_A.pdb", config["cap_flag"],
            max_cb_tilt_deg=config["max_cb_tilt_deg"]
        )
        anti_length, antiparallel_a = _build_aligned_peptide_template(
            sequence_a, [-139, 135], "anti_A.pdb", config["cap_flag"],
            max_cb_tilt_deg=None
        )
        if length != anti_length:
            raise ValueError(
                "Parallel and antiparallel templates have different lengths"
            )

        if not sequence_b:
            result = packpep_parallel_antiparallel(
                parallel_a,
                antiparallel_a,
                length,
                config["hybrid_type"],
                config["y_shift"],
                config["chains"],
                config["format_flag"],
                output_file,
                core_residues=config["core_residues"],
                strand_dist=config["x_spacing"],
                sheet_dist=config["z_spacing"],
                sheet_shift=config["x_shift"],
                registry_direction=config["registry_direction"],
            )
        else:
            para_length_b, parallel_b = _build_aligned_peptide_template(
                sequence_b, [-119, 113], "para_B.pdb", config["cap_flag"],
                max_cb_tilt_deg=config["max_cb_tilt_deg"]
            )
            anti_length_b, antiparallel_b = _build_aligned_peptide_template(
                sequence_b, [-139, 135], "anti_B.pdb", config["cap_flag"],
                max_cb_tilt_deg=None
            )
            template_lengths = {
                length,
                anti_length,
                para_length_b,
                anti_length_b,
            }
            if len(template_lengths) != 1:
                raise ValueError("The hybrid peptide templates have different lengths")
            result = packpep_parallel_antiparallel_two_component(
                parallel_a,
                parallel_b,
                antiparallel_a,
                antiparallel_b,
                length,
                config["hybrid_type"],
                config["y_shift"],
                pattern1,
                pattern2,
                config["format_flag"],
                output_file,
                core_residues=config["core_residues"],
                strand_dist=config["x_spacing"],
                sheet_dist=config["z_spacing"],
                sheet_shift=config["x_shift"],
                registry_direction=config["registry_direction"],
            )

    output_path = Path(output_file).resolve()
    if mode == "uniform":
        print(f"Sheets structure: Class {config['classes']}")
    else:
        print(f"Sheets structure: Hybrid {config['hybrid_type']}")

    if sequence_b:
        print(f'Sequence A: "{sequence_a}"')
        print(f'Sequence B: "{sequence_b}"')
        if len(pattern1) == len(pattern2):
            print(f"Strands per sheet: {len(pattern1)}")
        else:
            print(f"Strands in sheet 1: {len(pattern1)}")
            print(f"Strands in sheet 2: {len(pattern2)}")
    else:
        print(f'Sequence: "{sequence_a}"')
        print(f"Strands per sheet: {config['chains']}")

    print(f"Atoms: {len(result)}")
    print(f"Output PDB: {output_path}")
    return result


# %% Cell fcf611e5-e4a4-487c-914e-a4ce74087bdb
class RawDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter,
):
    pass


def read_arguments(
    argv: Optional[Sequence[str]] = None
) -> argparse.Namespace:
    """Function:
        Read inputs from the terminal.

    Parameters
    ----------
    argv : list of str or None, default=None
        Argument list used for testing. None reads the real terminal command line.

    Returns
    -------
    argparse.Namespace
        Parsed and conditionally validated sheet-builder inputs.
    """
    parser = argparse.ArgumentParser(
        description=("Initial Structure Builder v1.23: build beta sheets."),
        formatter_class=RawDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--seq",
        "--seqA",
        dest="sequence_a",
        required=True,
        default=argparse.SUPPRESS,
        help="Input peptide sequence. Or sequence-A when sequence-B is given.",
    )
    parser.add_argument(
        "-b",
        "--seqB",
        dest="sequence_b",
        default=argparse.SUPPRESS,
        help="Optional peptide sequence-B. It must have the same length as peptide A for this code verison.",
    )
    parser.add_argument(
        "-p1",
        "--pattern1",
        default=argparse.SUPPRESS,
        help="Sheet-1 A/B pattern. Required when peptide B is given.",
    )
    parser.add_argument(
        "-p2",
        "--pattern2",
        default=argparse.SUPPRESS,
        help="Sheet-2 A/B pattern. Required when peptide B is given.",
    )

    geometry = parser.add_mutually_exclusive_group(required=True)
    geometry.add_argument(
        "-c",
        "--class",
        dest="classes",
        type=int,
        choices=range(1, 11),
        default=argparse.SUPPRESS,
        help="Steric-zipper class from 1 to 10.",
    )
    geometry.add_argument(
        "-H",
        "--hybrid",
        dest="hybrid_type",
        type=int,
        choices=range(1, 7),
        default=argparse.SUPPRESS,
        help="Parallel-antiparallel hybrid sheet type from 1 to 6.",
    )

    parser.add_argument(
        "-dx",
        "--dx",
        type=float,
        default=4.8,
        help="Strand-strand distance along x, in Angstrom.",
    )
    parser.add_argument(
        "-dz",
        "--dz",
        type=float,
        default=11.5,
        help="Sheet-sheet distance along z, in Angstrom.",
    )
    parser.add_argument(
        "-x",
        "--x",
        type=float,
        default=0.0,
        help="Sheet-2 shifts in x direction in unit of [0.5*dx].",
    )
    parser.add_argument(
        "-y",
        "--y",
        type=float,
        default=0.0,
        help="Sheet-2 shifts in y direction in unit of residue spacing (3.465 Angstrom).",
    )
    parser.add_argument(
        "-n",
        "--chains",
        type=int,
        default=argparse.SUPPRESS,
        help="Strands per sheet for a one-sequence build; default is 8.",
    )
    parser.add_argument(
        "-r",
        "--core",
        choices=("e", "o"),
        default="e",
        help="packing e = even number packed inside, o = odd number packed inside.",
    )
    parser.add_argument(
        "-g",
        "--registry",
        choices=("minus", "plus"),
        default="minus",
        help=(
            "Relative direction of registry shift.\n"
            "Plus:\n"
            "Strand 1: 1 2 3 4 5 6\n"
            "Strand 2:   6 5 4 3 2 1\n"
            "Minus:\n"
            "Strand 1:   1 2 3 4 5 6\n"
            "Strand 2: 6 5 4 3 2 1\n"
        ),
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="format_flag",
        type=int,
        choices=(0, 1),
        default=0,
        help="PDB format: 0 = PepAD format; 1 = AMBER format.",
    )
    parser.add_argument(
        "-C",
        "--cap",
        dest="cap_flag",
        type=int,
        choices=(0, 1, 2),
        default=0,
        help="Terminal caps: 0 uncapped, 1 ACE+NME, or 2 ACE+NHE.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="peptide_sheets",
        help="Output PDB filename.",
    )
    parser.add_argument(
        "-t",
        "--tilt",
        dest="max_cb_tilt",
        type=float,
        default=10.0,
        help="Maximum angle between the middle-residue CA->CB vector and z-axis, in degrees.",
    )
    parser.add_argument(
        "--no-tilt",
        dest="max_cb_tilt",
        action="store_const",
        const=None,
        default=argparse.SUPPRESS,
        help="Unlimit the CA->CB vector and z-axis angle.",
    )

    args = parser.parse_args(argv)
    args.sequence_b = getattr(args, "sequence_b", None)
    args.pattern1 = getattr(args, "pattern1", None)
    args.pattern2 = getattr(args, "pattern2", None)
    args.classes = getattr(args, "classes", None)
    args.hybrid_type = getattr(args, "hybrid_type", None)
    chains_given = hasattr(args, "chains")
    has_sequence_b = args.sequence_b is not None
    has_pattern1 = args.pattern1 is not None
    has_pattern2 = args.pattern2 is not None

    if has_sequence_b and not (has_pattern1 and has_pattern2):
        parser.error("--seqB requires both --pattern1 and --pattern2")
    if not has_sequence_b and (has_pattern1 or has_pattern2):
        parser.error("--pattern1 and --pattern2 require --seqB")
    if not has_sequence_b:
        args.chains = args.chains if chains_given else 8
    elif chains_given:
        parser.error(
            "--chains is not used with --seqB; pattern lengths set strand counts"
        )
    return args


def main(argv=None):
    """Function:
        Build a peptide sheet from the concise terminal inputs.

    Parameters
    ----------
    argv : list of str or None, default=None
        Argument list used for testing. None reads the real terminal command line.

    Returns
    -------
    pandas.DataFrame
        Packed atoms from both sheets. The function also writes the output PDB file.
    """
    args = read_arguments(argv)
    return build_sheets(
        args.sequence_a,
        args.output,
        sequence_b=args.sequence_b,
        pattern1=args.pattern1,
        pattern2=args.pattern2,
        classes=args.classes,
        hybrid_type=args.hybrid_type,
        core_residues=args.core,
        chains=getattr(args, "chains", 8),
        format_flag=args.format_flag,
        cap_flag=args.cap_flag,
        y_shift=args.y,
        registry_direction=args.registry,
        x_spacing=args.dx,
        z_spacing=args.dz,
        x_shift=args.x,
        max_cb_tilt_deg=args.max_cb_tilt,
    )


def cli() -> None:
    """Run the installed ``builder`` terminal command."""
    main()


if __name__ == "__main__" and "get_ipython" not in globals():
    cli()
