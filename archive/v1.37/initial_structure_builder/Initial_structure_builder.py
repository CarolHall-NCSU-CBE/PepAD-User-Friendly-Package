#!/usr/bin/env python
# coding: utf-8

# In[58]:


import numpy as np
# from PeptideBuilder import Geometry
# import PeptideBuilder
import pandas as pd
import scipy.linalg
import matplotlib.pyplot as plt
import argparse
# import Bio.PDB
# import hydride
# import subprocess


# In[59]:


def linear_fitting(data):
#      reference:    
    # https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d

    datamean = data.mean(axis=0)
    uu, dd, vv = np.linalg.svd(data - datamean)
    return vv[0]

def linear_fitting_3D_points(points):

#      reference:
#      https://www.doc88.com/p-8189740853644.html  
#     《三维空间点中基于最小二乘法的分段直线拟合方法》 薛丽红，2015年7月，齐齐哈尔学报，第31卷第4期

    Sum_X=0.0
    Sum_Y=0.0
    Sum_Z=0.0
    Sum_XZ=0.0
    Sum_YZ=0.0
    Sum_Z2=0.0
    n=len(points)
    for i in range(0,len(points)):
        xi=points[i][0]
        yi=points[i][1]
        zi=points[i][2]

        Sum_X = Sum_X + xi
        Sum_Y = Sum_Y + yi
        Sum_Z = Sum_Z + zi
        Sum_XZ = Sum_XZ + xi*zi
        Sum_YZ = Sum_YZ + yi*zi
        Sum_Z2 = Sum_Z2 + zi**2

    den = n*Sum_Z2 - Sum_Z * Sum_Z # 公式分母
    k1 = (n*Sum_XZ - Sum_X * Sum_Z)/ den
    b1 = (Sum_X - k1 * Sum_Z)/n
    k2 = (n*Sum_YZ - Sum_Y * Sum_Z)/ den
    b2 = (Sum_Y - k2 * Sum_Z)/n
    
    return k1, b1, k2, b2

def unit_vector(vector, label='vector', allow_zero=False):
    """Return a normalized vector; optionally return None for near-zero input."""
    vector = np.asarray(vector, dtype=np.float64)
    norm = np.linalg.norm(vector)
    if norm < 1e-8:
        if allow_zero:
            return None
        raise ValueError(f'Cannot normalize near-zero {label}.')
    return vector / norm
    
def cross2d(a, b):
    return a[0]*b[1] - a[1]*b[0]
    
def rotation_angle(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    
    cos_angle = np.dot(v1_u, v2_u)/(np.sqrt(np.dot(v1_u, v1_u))*np.sqrt(np.dot(v2_u, v2_u)))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    angle = np.arccos(cos_angle)
    
    # judge rotation direction
    cross_product = cross2d(v1_u, v2_u) # cross_product OF 2d vector, it return a scalar
    if cross_product < 0:
        direction = -1
    else:
        direction = 1
    
    angle = angle * direction
    return angle

def rotation_coordinates(df, axis, angle):
    df_temp = df.copy()
    data = df_temp[['x', 'y', 'z']].to_numpy()

    if(axis == 'x'):
        M = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
    elif (axis == 'y'):
        M = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
    elif (axis == 'z'):
        M = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    elif (axis == "none"):
        return df_temp
    
    rotated_data = (M @ data.T).T  # Transpose to apply matrix multiplication correctly
    df_temp[['x', 'y', 'z']] = rotated_data

    return df_temp

def add_NH_remove_OH(pdbfile, o_filename):
    pdb=[]
    columns = ['atom', 'anum', 'atom_name', 'aa_name', 'resid', 'x', 'y', 'z']

    with open(pdbfile, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("ATOM"):
            components = line.split()

            atom = line[0:6].strip()           # ATOM
            anum = int(line[6:11])            # Atom number
            atom_name = line[11:17].strip()   # Atom name
            aa_name = line[17:21].strip()     # Amino acid name
            resid = int(line[22:26])          # Residue ID
            x = float(line[26:38])            # X-coordinate
            y = float(line[38:46])            # Y-coordinate
            z = float(line[46:54])            # Z-coordinate

            # Create a dictionary for the current line
            row = {
                'atom': atom,
                'anum': anum,
                'atom_name': atom_name,
                'aa_name': aa_name,
                'resid': resid,
                'x': x,
                'y': y,
                'z': z
            }

            # Append the dictionary to the list
            pdb.append(row)
            
    df = pd.DataFrame(pdb, columns=columns)
    ##### remove "H" on "OXT" ############### 
    df = df[df['atom_name'] != "HXT"]
    df.loc[(df['atom_name'] == 'H') & (df['resid'] == 1), 'atom_name'] = 'H1'
    
    h2 = df[df['atom_name'] == 'H2'].iloc[0].copy()       
    h3 = h2.copy()
    h3['atom_name'] = 'H3'
    idx = df[df['atom_name'] == 'H2'].index[0] + 1                        # index of H2 in the df
    df = pd.concat([df.iloc[:idx], pd.DataFrame([h3]), df.iloc[idx:]]).reset_index(drop=True) # combine dataframe
    df['anum'] = range(1, len(df) + 1)                                                        # Reindex atom number

    ################## calculate H coordinates #########
    # Solve the corrected system of equations
    if df[df['resid'] == 1].iloc[0]['aa_name'] == 'GLY':
        cb = df[(df['atom_name'] == 'HA2') & (df['resid'] == 1)].iloc[0][['x', 'y', 'z']].to_numpy()
    else:
        cb = df[(df['atom_name'] == 'CB') & (df['resid'] == 1)].iloc[0][['x', 'y', 'z']].to_numpy()

    ca = df[(df['atom_name'] == 'CA') & (df['resid'] == 1)].iloc[0][['x', 'y', 'z']].to_numpy()
    n = df[(df['atom_name'] == 'N') & (df['resid'] == 1)].iloc[0][['x', 'y', 'z']].to_numpy()

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

    df.loc[df['atom_name'] == 'H1', ['x', 'y', 'z']] = h_coordinate
    df.loc[df['atom_name'] == 'H2', ['x', 'y', 'z']] = rotated_h_coords[0]
    df.loc[df['atom_name'] == 'H3', ['x', 'y', 'z']] = rotated_h_coords[1]

    with open(o_filename, 'w') as f:
        for index, row in df.iterrows():
            if row['resid'] == 1:
                row['aa_name'] = "N" + row['aa_name']
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])
                
            elif row['resid'] == df['resid'].max():
                row['aa_name'] = "C" + row['aa_name']
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])
                
            else:
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])
            
            f.write(output)

def add_caps(pdbfile, o_filename, type_flag):            
    # C-NME bond length = 1.291


    NHE_data = [
        {"anum": 1, "atom_name": "N",    "aa_name": "NHE", "resid": 1, "x": 2.194, "y": 1.598, "z": -0.000},
        {"anum": 2, "atom_name": "HN1",  "aa_name": "NHE", "resid": 1, "x": 3.035, "y": 1.039, "z": -0.000},
        {"anum": 3, "atom_name": "HN2",  "aa_name": "NHE", "resid": 1, "x": 2.250, "y": 2.606, "z": -0.000},
    ]
    
    NME_data = [
        {"anum": 1, "atom_name": "N",    "aa_name": "NME", "resid": 1, "x": 3.326, "y": 1.548, "z": -0.000},
        {"anum": 2, "atom_name": "H",    "aa_name": "NME", "resid": 1, "x": 3.909, "y": 0.724, "z": -0.000},
        {"anum": 3, "atom_name": "CH3",  "aa_name": "NME", "resid": 1, "x": 3.970, "y": 2.846, "z": -0.000},
        {"anum": 4, "atom_name": "HH31", "aa_name": "NME", "resid": 1, "x": 3.212, "y": 3.629, "z":  0.000},
        {"anum": 5, "atom_name": "HH32", "aa_name": "NME", "resid": 1, "x": 4.592, "y": 2.943, "z":  0.890},
        {"anum": 6, "atom_name": "HH33", "aa_name": "NME", "resid": 1, "x": 4.592, "y": 2.943, "z": -0.890},
    ]
    
    ACE_data = [
        {"anum": 1, "atom_name": "H1",   "aa_name": "ACE", "resid": 1, "x": 2.000, "y": 1.000, "z": -0.000},
        {"anum": 2, "atom_name": "CH3",  "aa_name": "ACE", "resid": 1, "x": 2.000, "y": 2.090, "z":  0.000},
        {"anum": 3, "atom_name": "H2",   "aa_name": "ACE", "resid": 1, "x": 1.486, "y": 2.454, "z":  0.890},
        {"anum": 4, "atom_name": "H3",   "aa_name": "ACE", "resid": 1, "x": 1.486, "y": 2.454, "z": -0.890},
        {"anum": 5, "atom_name": "C",    "aa_name": "ACE", "resid": 1, "x": 3.427, "y": 2.641, "z": -0.000},
        {"anum": 6, "atom_name": "O",    "aa_name": "ACE", "resid": 1, "x": 4.391, "y": 1.877, "z": -0.000},
    ] 

    NHE = pd.DataFrame(NHE_data)
    NME = pd.DataFrame(NME_data)
    ACE = pd.DataFrame(ACE_data)
    
    pdb=[]
    columns = ['atom', 'anum', 'atom_name', 'aa_name', 'resid', 'x', 'y', 'z']

    with open(pdbfile, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("ATOM"):
            components = line.split()

            atom = line[0:6].strip()          # ATOM
            anum = int(line[6:11])            # Atom number
            atom_name = line[11:17].strip()   # Atom name
            aa_name = line[17:21].strip()     # Amino acid name
            resid = int(line[22:26])          # Residue ID
            x = float(line[26:38])            # X-coordinate
            y = float(line[38:46])            # Y-coordinate
            z = float(line[46:54])            # Z-coordinate

            row = {
                'atom': atom,
                'anum': anum,
                'atom_name': atom_name,
                'aa_name': aa_name,
                'resid': resid,
                'x': x,
                'y': y,
                'z': z
            }

            pdb.append(row)
    
    df = pd.DataFrame(pdb, columns=columns)
    
##### remove "H" on "OXT" ############### 
    
    
    df = df[(df['atom_name'] != "HXT")]                      # remove H on OXT

    # Keep one backbone N-H for the ACE amide bond and normalize its name to "H".
    # Hydride/PDB writers may call these H, H1/H2, or HN/HN1/HN2; Gly at the
    # N-terminus is fine as long as one N-terminal backbone H is present.
    n_term_h_names = ['H', 'H1', 'HN', 'HN1', 'H2', 'HN2', 'H3', 'HN3']
    n_term_h = df[(df['resid'] == 1) & (df['atom_name'].isin(n_term_h_names))]
    if n_term_h.empty:
        raise ValueError("Cannot add ACE cap: no N-terminal backbone H was found on residue 1.")
    keep_h_idx = n_term_h.index[0]
    remove_h_idx = n_term_h.index[1:]
    df = df.drop(remove_h_idx)
    df.loc[keep_h_idx, 'atom_name'] = 'H'

    df['resid'] += 1                                         # increase residue ID by 1
    
    ####### add ACE group ###########
    df = pd.concat([ACE, df], ignore_index=True)

    ###### calcualte NME or NHE position ########
    last_resid = df[df['atom_name'] == 'C']['resid'].max()
    C = df[(df['atom_name'] == 'C') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
    OXT = df[(df['atom_name'] == 'OXT') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
    C_OXT_vector = OXT - C
    C_N_vector = C_OXT_vector * 1.291 / np.linalg.norm(C_OXT_vector) 
    N = C + C_N_vector
    df = df[~((df['atom_name'] == "OXT") & (df['resid'] == last_resid))]          # remove OXT
    
    if type_flag == 1:
        n_pos = NME[NME['atom_name'] == "N"][['x', 'y', 'z']].values[0]
        ####### add NME group ###########
        NME['x'] = NME['x'] - n_pos[0] + N[0]
        NME['y'] = NME['y'] - n_pos[1] + N[1]
        NME['z'] = NME['z'] - n_pos[2] + N[2]
        NME['resid'] += last_resid
    
        df = pd.concat([df, NME], ignore_index=True)
        df['anum'] = range(1, len(df) + 1) 
        
        ####### Rotate NME group #######
        CA= df[(df['atom_name'] == 'CA') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
        N = df[(df['aa_name'] == 'NME') & (df['atom_name'] == 'N')].iloc[0][['x', 'y', 'z']].to_numpy()
        H = df[(df['aa_name'] == 'NME') & (df['atom_name'] == 'H')].iloc[0][['x', 'y', 'z']].to_numpy()     # H from NME
        CH3 = df[(df['aa_name'] == 'NME') & (df['atom_name'] == 'CH3')].iloc[0][['x', 'y', 'z']].to_numpy()
        O = df[(df['resid'] == last_resid) & (df['atom_name'] == 'O')].iloc[0][['x', 'y', 'z']].to_numpy()
        C = np.asarray(C, dtype=np.float64); CA= np.asarray(CA, dtype=np.float64); N = np.asarray(N, dtype=np.float64); 
        H = np.asarray(H, dtype=np.float64); CH3 = np.asarray(CH3, dtype=np.float64); O = np.asarray(O, dtype=np.float64);
        
        v1 = C-N; v2 = CH3 - N
        v3 = rotate_vector_plane_3D(v1, v2, np.pi*2/3) * np.linalg.norm(v2)      # rotate v2 to v1 by 120 degree
        v4 = rotate_vector_plane_3D(v1, v2, np.pi*4/3) * np.linalg.norm(H - N)      # also get position of H
    
        CH3 = N + v3                                         # correct position of CH3, N
        H = N + v4  
        
        df.loc[(df['atom_name'] == 'CH3') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = CH3
        df.loc[(df['atom_name'] == 'H') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = H
        
        ################## calculate H coordinates #########
        p1 = H
        p2 = N
        p3 = CH3
    
        # Parameters
        dihedral_angle_deg = 97.150; bond_angle_deg = 110.00; bond_length = 1.078
    
        h_coordinate = place_fourth_atom(p1, p2, p3, bond_length, bond_angle_deg, dihedral_angle_deg) ## important
        
        angles = [120, 240]  # Angles to rotate by
        rotated_h_coords = rotate_around_axis(p2, p3, h_coordinate, angles)
        
        df.loc[(df['atom_name'] == 'HH31') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = h_coordinate
        df.loc[(df['atom_name'] == 'HH32') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = rotated_h_coords[0]
        df.loc[(df['atom_name'] == 'HH33') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = rotated_h_coords[1]
        
        ##################### rotate C(amino acid) - N(NME) axis to 180 degree ################  
        theta = dihedral_angle(O, C, N, CH3)
        rotation_angle = 0 * np.pi - theta
    
        # Define rotation axis (C → N)
        rotation_axis = N-C
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize axis
        rotation_center = N
    
        # Copy NME group and rotate its atoms
        NME_atoms = df[df['aa_name'] == "NME"].copy()
    
        for i, row in NME_atoms.iterrows():
            if row['atom_name'] != "N":
                v = np.array(row[['x', 'y', 'z']].values, dtype=np.float64) - rotation_center  # Translate to origin
                v_rot = rotate_vector(v, rotation_axis, rotation_angle)  # Rotate
                NME_atoms.loc[i, ['x', 'y', 'z']] = (v_rot + rotation_center).tolist()  # Translate back
    
        # Replace old NME with rotated version
        df = df[df['aa_name'] != "NME"]  
        df = pd.concat([df, NME_atoms], ignore_index=True)

    
    ##################### rotate CA(amino acid) - C(Amino acid) axis to 180 degree ################

        terminal_ha = df[(df['resid'] == last_resid) & (df['atom_name'].isin(['HA', 'HA2', 'HA3']))]
        if terminal_ha.empty:
            raise ValueError(f"Cannot orient C-terminal cap: no HA/HA2/HA3 found on residue {last_resid}.")
        HA= terminal_ha.iloc[0][['x', 'y', 'z']].to_numpy()
        CA= df[(df['atom_name'] == 'CA') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
        HA = np.asarray(HA, dtype=np.float64); CA= np.asarray(CA, dtype=np.float64)
    
        theta = dihedral_angle(HA, CA, C, O)
        rotation_angle = (-164.77)/180 * np.pi - theta
    
        # Define rotation axis
        rotation_axis = C-CA
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize axis
        rotation_center = C
    
        # Copy NME group and rotate its atoms
        rotating_atoms = df[(df['aa_name'] == "NME") | ((df['atom_name'] == 'O') & (df['resid'] == last_resid))].copy()
    
    
        for i, row in rotating_atoms.iterrows():
            v = np.array(row[['x', 'y', 'z']].values, dtype=np.float64) - rotation_center  # Translate to origin
            v_rot = rotate_vector(v, rotation_axis, rotation_angle)  # Rotate
            rotating_atoms.loc[i, ['x', 'y', 'z']] = (v_rot + rotation_center).tolist()  # Translate back
    
        # Replace old NME with rotated version
        df = df[~((df['aa_name'] == "NME") | ((df['atom_name'] == 'O') & (df['resid'] == last_resid)))] 
        df = pd.concat([df, rotating_atoms], ignore_index=True)

    elif type_flag == 2:
        n_pos = NHE[NHE['atom_name'] == "N"][['x', 'y', 'z']].values[0]
        ####### add NHE group ###########
        NHE['x'] = NHE['x'] - n_pos[0] + N[0]
        NHE['y'] = NHE['y'] - n_pos[1] + N[1]
        NHE['z'] = NHE['z'] - n_pos[2] + N[2]
        NHE['resid'] += last_resid
    
        df = pd.concat([df, NHE], ignore_index=True)
        df['anum'] = range(1, len(df) + 1) 
        ####### Rotate NME group #######
        CA= df[(df['atom_name'] == 'CA') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
        HN1 = df[(df['aa_name'] == 'NHE') & (df['atom_name'] == 'HN1')].iloc[0][['x', 'y', 'z']].to_numpy()
        HN2 = df[(df['aa_name'] == 'NHE') & (df['atom_name'] == 'HN2')].iloc[0][['x', 'y', 'z']].to_numpy()    
        O = df[(df['resid'] == last_resid) & (df['atom_name'] == 'O')].iloc[0][['x', 'y', 'z']].to_numpy()
        
        C = np.asarray(C, dtype=np.float64); CA= np.asarray(CA, dtype=np.float64); N = np.asarray(N, dtype=np.float64); 
        HN1 = np.asarray(HN1, dtype=np.float64); HN2 = np.asarray(HN2, dtype=np.float64); O = np.asarray(O, dtype=np.float64);
        
        v1 = C-N; v2 = HN1 - N
        v3 = rotate_vector_plane_3D(v1, v2, np.pi*2/3) * np.linalg.norm(v2)           # rotate v2 to v1 by 120 degree
        v4 = rotate_vector_plane_3D(v1, v2, np.pi*4/3) * np.linalg.norm(HN2 - N)      # also get position of H
    
        HN1 = N + v3                                         # correct position of CH3, N
        HN2 = N + v4  
        
        df.loc[(df['atom_name'] == 'HN1') & (df['aa_name'] == 'NHE'), ['x', 'y', 'z']] = HN1
        df.loc[(df['atom_name'] == 'HN2') & (df['aa_name'] == 'NHE'), ['x', 'y', 'z']] = HN2
        
        # ################## calculate H coordinates #########
        # p1 = H
        # p2 = N
        # p3 = CH3
    
        # # Parameters
        # dihedral_angle_deg = 97.150; bond_angle_deg = 110.00; bond_length = 1.078
    
        # h_coordinate = place_fourth_atom(p1, p2, p3, bond_length, bond_angle_deg, dihedral_angle_deg) ## important
        
        # angles = [120, 240]  # Angles to rotate by
        # rotated_h_coords = rotate_around_axis(p2, p3, h_coordinate, angles)
        
        # df.loc[(df['atom_name'] == 'HH31') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = h_coordinate
        # df.loc[(df['atom_name'] == 'HH32') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = rotated_h_coords[0]
        # df.loc[(df['atom_name'] == 'HH33') & (df['aa_name'] == 'NME'), ['x', 'y', 'z']] = rotated_h_coords[1]
        
        ##################### rotate C(amino acid) - N(NME) axis to 180 degree ################  
        theta = dihedral_angle(O, C, N, HN1)
        rotation_angle = 0 * np.pi - theta
    
        # Define rotation axis (C → N)
        rotation_axis = N-C
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize axis
        rotation_center = N
    
        # Copy NME group and rotate its atoms
        NHE_atoms = df[df['aa_name'] == "NHE"].copy()
    
        for i, row in NHE_atoms.iterrows():
            if row['atom_name'] != "N":
                v = np.array(row[['x', 'y', 'z']].values, dtype=np.float64) - rotation_center  # Translate to origin
                v_rot = rotate_vector(v, rotation_axis, rotation_angle)  # Rotate
                NHE_atoms.loc[i, ['x', 'y', 'z']] = (v_rot + rotation_center).tolist()  # Translate back
    
        # Replace old NME with rotated version
        df = df[df['aa_name'] != "NHE"]  
        df = pd.concat([df, NHE_atoms], ignore_index=True)

    
    ##################### rotate CA(amino acid) - C(Amino acid) axis to 180 degree ################

        terminal_ha = df[(df['resid'] == last_resid) & (df['atom_name'].isin(['HA', 'HA2', 'HA3']))]
        if terminal_ha.empty:
            raise ValueError(f"Cannot orient C-terminal cap: no HA/HA2/HA3 found on residue {last_resid}.")
        HA= terminal_ha.iloc[0][['x', 'y', 'z']].to_numpy()
        CA= df[(df['atom_name'] == 'CA') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
        HA = np.asarray(HA, dtype=np.float64); CA= np.asarray(CA, dtype=np.float64)
    
        theta = dihedral_angle(HA, CA, C, O)
        rotation_angle = (-164.77)/180 * np.pi - theta
    
        # Define rotation axis
        rotation_axis = C-CA
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize axis
        rotation_center = C
    
        # Copy NME group and rotate its atoms
        rotating_atoms = df[(df['aa_name'] == "NHE") | ((df['atom_name'] == 'O') & (df['resid'] == last_resid))].copy()
    
    
        for i, row in rotating_atoms.iterrows():
            v = np.array(row[['x', 'y', 'z']].values, dtype=np.float64) - rotation_center  # Translate to origin
            v_rot = rotate_vector(v, rotation_axis, rotation_angle)  # Rotate
            rotating_atoms.loc[i, ['x', 'y', 'z']] = (v_rot + rotation_center).tolist()  # Translate back
    
        # Replace old NME with rotated version
        df = df[~((df['aa_name'] == "NHE") | ((df['atom_name'] == 'O') & (df['resid'] == last_resid)))] 
        df = pd.concat([df, rotating_atoms], ignore_index=True)
    
    ###############################################################################################################################


    
    ############## correct ACE position ##############
    CA= df[(df['atom_name'] == 'CA') & (df['resid'] == 2)].iloc[0][['x', 'y', 'z']].to_numpy()
    N = df[(df['atom_name'] == 'N') & (df['resid'] == 2)].iloc[0][['x', 'y', 'z']].to_numpy()
    H = df[(df['atom_name'] == 'H') & (df['resid'] == 2)].iloc[0][['x', 'y', 'z']].to_numpy()
    CA= np.asarray(CA, dtype=np.float64); N = np.asarray(N, dtype=np.float64); H = np.asarray(H, dtype=np.float64);
    
    v1 = H-N; v2 = CA - N
    v3 = rotate_vector_plane_3D(v1, v2, np.pi*2/3) * 1.290      # rotate v2 to v1 by 120 degree
    C = N + v3                                                  # correct position of CH3, N
  
    df.loc[(df['atom_name'] == 'C') & (df['aa_name'] == 'ACE'), ['x', 'y', 'z']] = C
    
    v1 = N - C
    # Use the CA-N direction as the ACE reference; H-N swaps the ACE O/CH3 sides.
    v2 = CA - N
    v3 = rotate_vector_plane_3D(v1, v2, np.pi*4/3) * 1.238
    v4 = rotate_vector_plane_3D(v1, v2, np.pi*2/3) * 1.514
    O = C + v3
    CH3 = C + v4
    df.loc[(df['atom_name'] == 'O') & (df['aa_name'] == 'ACE'), ['x', 'y', 'z']] = O
    df.loc[(df['atom_name'] == 'CH3') & (df['aa_name'] == 'ACE'), ['x', 'y', 'z']] = CH3
    
    ####### place H3 #########
    p1 = O
    p2 = C
    p3 = CH3

    # Parameters
    dihedral_angle_deg = 60
    bond_angle_deg = 110.00
    bond_length = 1.080

    h_coordinate = place_fourth_atom(p1, p2, p3, bond_length, bond_angle_deg, dihedral_angle_deg)
    
    angles = [120, 240]  # Angles to rotate by
    rotated_h_coords = rotate_around_axis(p2, p3, h_coordinate, angles)
    
    df.loc[(df['atom_name'] == 'H1') & (df['aa_name'] == 'ACE'), ['x', 'y', 'z']] = h_coordinate
    df.loc[(df['atom_name'] == 'H2') & (df['aa_name'] == 'ACE'), ['x', 'y', 'z']] = rotated_h_coords[0]
    df.loc[(df['atom_name'] == 'H3') & (df['aa_name'] == 'ACE'), ['x', 'y', 'z']] = rotated_h_coords[1]

    ####### output structure ###########
    with open(o_filename, 'w') as f:
        for index, row in df.iterrows():
            if len(row['aa_name'])==4:
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])
            else:
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])

            f.write(output)



# In[60]:


def rotate_vector(v, axis, angle):
    """
    Rotate vector `v` around `axis` by `angle` radians using Rodrigues' rotation formula.
    """
    axis = axis / np.linalg.norm(axis)  # Normalize the axis
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    cross_prod = np.cross(axis, v)
    dot_prod = np.dot(axis, v)
    return v * cos_theta + cross_prod * sin_theta + axis * dot_prod * (1 - cos_theta)

def rotate_around_axis(p1, p2, p3, angles):
    """
    Rotate H around the p1->p2 axis by given angles.
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
    """
    Rotate multiple points around axis defined by (p1 -> p2) by `angle` radians.
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
    
    n = np.cross(u, v)
    n = np.asarray(n, dtype=np.float64); 

    n /= np.linalg.norm(n)  # Normalize
    v_rot = v * np.cos(angle) + np.cross(n, v) * np.sin(angle) + n * np.dot(n, v) * (1 - np.cos(angle))
    v_rot /= np.linalg.norm(v_rot)
    
    return v_rot


def dihedral_angle(p1, p2, p3, p4):
    """Calculate the dihedral angle between four points.

    Args:
        p1, p2, p3, p4 (numpy.ndarray): Coordinates of the four points.

    Returns:
        float: The dihedral angle in degrees.
    """

    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    cos_phi = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
    sin_phi = np.dot(np.cross(n1, n2), b2) / (np.linalg.norm(b2) * np.linalg.norm(n1) * np.linalg.norm(n2))

    phi = np.arctan2(sin_phi, cos_phi)

    return phi


def place_fourth_atom(p1, p2, p3, bond_length, bond_angle_deg, dihedral_deg):
    """
    p1 = CB coordinates (x1, y1, z1)
    p2 = CA coordinates (x2, y2, z2)
    p3 = N coordinates (x3, y3, z3)
    p4 = H (x, y, z) unknown
    bond_length = length of N-H bond
    bond_angle_deg = CA-N-H angle (in degrees)
    dihedral_deg = CB-CA-N-H dihedral angle (in degrees)
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
    y = -bond_length * np.sin(bond_angle) * np.cos(dihedral + np.pi/2)
    z = bond_length * np.sin(bond_angle) * np.sin(dihedral + np.pi/2)

    # Position of H in global coordinates
    p4 = p3 + x * e1 + y * e2 + z * e3
    return p4


# In[61]:


def _read_simple_pdb(pdbfile):
    records = []
    with open(pdbfile, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            records.append({
                'atom': line[0:6].strip() or 'ATOM',
                'anum': int(line[6:11]),
                'atom_name': line[11:17].strip(),
                'aa_name': line[17:21].strip(),
                'resid': int(line[22:26]),
                'x': float(line[26:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
            })
    return pd.DataFrame(records, columns=['atom', 'anum', 'atom_name', 'aa_name', 'resid', 'x', 'y', 'z'])


def _write_simple_pdb(df, o_filename):
    with open(o_filename, 'w') as f:
        for _, row in df.iterrows():
            f.write("ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                int(row['anum']), row['atom_name'], row['aa_name'], int(row['resid']),
                float(row['x']), float(row['y']), float(row['z'])
            ))


def copy_pdb_file(input_file, output_file):
    with open(input_file, 'r') as src, open(output_file, 'w') as dst:
        dst.write(src.read())


def _place_opposite_bisector(center, neighbors, bond_length):
    direction = np.zeros(3, dtype=np.float64)
    for neighbor in neighbors:
        unit = unit_vector(neighbor - center, allow_zero=True)
        if unit is not None:
            direction += unit
    direction = unit_vector(-direction, allow_zero=True)
    if direction is None:
        return None
    return center + bond_length * direction





# In[62]:


# Local PeptideBuilder compatibility layer, vendored from PeptideBuilder 1.1.0.
# Original authors: Matthew Z. Tien, Dariya K. Sydykova, Austin G. Meyer,
# and Claus O. Wilke. The upstream files state that they are provided under
# the MIT License. This cell removes the external PeptideBuilder package
# dependency while preserving the same peptide-building geometry.

"""This module is part of the PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

The Geometry module contains the default geometries of
all 20 amino acids. The main function to be used is the
geometry() function, which returns the default geometry
for the requested amino acid.

This file is provided to you under the MIT License."""

import random
from typing import List


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
    """Generates the geometry of the requested amino acid.
    The amino acid needs to be specified by its single-letter
    code. If an invalid code is specified, the function
    returns the geometry of Glycine."""
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
import math, warnings
from typing import Optional, Union
import numpy as np


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
    if isinstance(value, Vector):
        return value.array
    if isinstance(value, Atom):
        return value.coord
    return np.asarray(value, dtype=float)


def rotaxis(theta, vector):
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
    a = _as_array(v1) - _as_array(v2)
    b = _as_array(v3) - _as_array(v2)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        raise ValueError("Cannot calculate angle with a zero-length vector.")
    return math.acos(np.clip(np.dot(a, b) / denom, -1.0, 1.0))


def calc_dihedral(v1, v2, v3, v4):
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
    """Creates a Glycine residue"""
    res = Residue((" ", segID, " "), "GLY", "    ")

    res.add(N)
    res.add(CA)
    res.add(C)
    res.add(O)
    return res


def makeAla(segID: int, N, CA, C, O, geo: AlaGeo) -> Residue:
    """Creates an Alanine residue"""
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
    """Creates a Serine residue"""
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
    """Creates a Cysteine residue"""
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
    """Creates a Valine residue"""
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
    """Creates an Isoleucine residue"""
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
    """Creates a Leucine residue"""
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
    """Creates a Threonine residue"""
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
    """Creates an Arginie residue"""
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
    """Creates a Lysine residue"""
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
    """Creates an Aspartic Acid residue"""
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
    """Creates an Asparagine residue"""
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
    """Creates a Glutamic Acid residue"""
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
    """Creates a Glutamine residue"""
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
    """Creates a Methionine residue"""
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
    """Creates a Histidine residue"""
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
    """Creates a Proline residue"""
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
    """Creates a Phenylalanine residue"""
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
    """Creates a Tyrosine residue"""
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
    """Creates a Tryptophan residue"""
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
    """Creates a new structure containing a single amino acid. The type and
    geometry of the amino acid are determined by the argument, which has to be
    either a geometry object or a single-letter amino acid code.
    The amino acid will be placed into chain A of model 0."""

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
    """Returns the last residue of chain A model 0 of the given structure.

    This function is a helper function that should not normally be called
    directly."""

    # If the following line doesn't work we're in trouble.
    # Likely initialize_res() wasn't called.
    resRef = structure[0]["A"].child_list[-1]

    # If the residue is not an amino acid we're in trouble.
    # Likely somebody is trying to append residues to an existing
    # structure that has non-amino-acid molecules in the chain.
    assert is_aa(resRef)

    return resRef


def add_residue_from_geo(structure: Structure, geo: Geo) -> Structure:
    """Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added is determined by
    the geometry object given as second argument.

    This function is a helper function and should not normally be called
    directly. Call add_residue() instead."""
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
    """Place a sequence of amino acids into a peptide in the extended
    conformation. The argument AA_chain holds the sequence of amino
    acids to be used."""
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
    """Adds a residue to chain A model 0 of the given structure, and
    returns the new structure. The residue to be added can be specified
    in two ways: either as a geometry object (in which case
    the remaining arguments phi, psi_im1, and omega are ignored) or as a
    single-letter amino-acid code. In the latter case, the optional
    arguments phi, psi_im1, and omega specify the corresponding backbone
    angles.

    When omega is specified, it needs to be a value greater than or equal
    to -360. Values below -360 are ignored."""

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
    """Place a sequence of amino acids into a peptide with specified
    backbone dihedral angles. The argument AA_chain holds the
    sequence of amino acids to be used. The arguments phi and psi_im1 hold
    lists of backbone angles, one for each amino acid, *starting from
    the second amino acid in the chain*. The argument
    omega (optional) holds a list of omega angles, also starting from
    the second amino acid in the chain."""
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
    """Creates a structure out of a list of geometry objects."""
    model_structure = initialize_res(geos[0])
    for i in range(1, len(geos)):
        add_residue(model_structure, geos[i])

    return model_structure


def add_terminal_OXT(structure: Structure, C_OXT_length: float = 1.23) -> Structure:
    """Adds a terminal oxygen atom ('OXT') to the last residue of chain A model 0 of the given structure, and returns the new structure. The OXT atom object will be contained in the last residue object of the structure.

This function should be used only when the structure object is completed and no further residues need to be appended."""

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




# In[63]:


from pathlib import Path

AA1_TO_ROTAMER = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIE', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}


def structure_to_simple_dataframe(structure):
    rows = []
    anum = 1
    for model in structure.child_list:
        for chain in model.child_list:
            for residue in chain.child_list:
                resid = residue.get_id()[1]
                for atom in residue.child_list:
                    x, y, z = atom.coord
                    rows.append({
                        'atom': 'ATOM',
                        'anum': anum,
                        'atom_name': atom.name,
                        'aa_name': residue.resname,
                        'resid': resid,
                        'x': float(x),
                        'y': float(y),
                        'z': float(z),
                    })
                    anum += 1
    return pd.DataFrame(rows, columns=['atom', 'anum', 'atom_name', 'aa_name', 'resid', 'x', 'y', 'z'])


def read_rotamer_template(residue_name, rotamer_dir='RotamerLibrary'):
    template_file = Path(rotamer_dir) / residue_name
    if not template_file.exists():
        raise FileNotFoundError(f"Rotamer template not found: {template_file}")
    return _read_simple_pdb(str(template_file))


def _kabsch_transform(moving_points, target_points):
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
    coords = np.asarray(coords, dtype=np.float64)
    return coords @ rotation.T + translation


def build_target_beta_backbone(sequence, angles):
    geos = []
    for i, aa in enumerate(sequence):
        geo = Geometry.geometry(aa)
        if i > 0:
            geo.phi = angles[0]
            geo.psi_im1 = angles[1]
        geos.append(geo)
    structure = PeptideBuilder.make_structure_from_geos(geos)
    PeptideBuilder.add_terminal_OXT(structure)
    return structure_to_simple_dataframe(structure)





# In[64]:


# Clean peptide builder: poly-A/G backbone + rotamer side-chain transplant

BACKBONE_TEMPLATE_ATOMS = {'N', 'H', 'H1', 'H2', 'H3', 'CA', 'HA', 'HA2', 'HA3', 'C', 'O', 'OXT'}


def _atom_xyz_from_residue(residue_df, atom_name):
    hit = residue_df[residue_df['atom_name'] == atom_name]
    if hit.empty:
        return None
    return hit.iloc[0][['x', 'y', 'z']].to_numpy(dtype=np.float64)


def _append_atom(rows, atom_name, aa_name, resid, xyz):
    rows.append({
        'atom': 'ATOM',
        'anum': 0,
        'atom_name': atom_name,
        'aa_name': aa_name,
        'resid': int(resid),
        'x': float(xyz[0]),
        'y': float(xyz[1]),
        'z': float(xyz[2]),
    })


def Add_backbone_hydrogen(backbone_df, sequence):
    # Add only peptide backbone and alpha hydrogens to an Ala/Gly scaffold.
    df = backbone_df.copy()
    sequence = sequence.upper().strip()
    additions = []
    residues = sorted(df['resid'].unique())
    residue_lookup = {resid: df[df['resid'] == resid].copy() for resid in residues}

    for i, resid in enumerate(residues):
        res = residue_lookup[resid]
        original_aa = sequence[i]
        aa_name = res.iloc[0]['aa_name']
        n = _atom_xyz_from_residue(res, 'N')
        ca = _atom_xyz_from_residue(res, 'CA')
        c = _atom_xyz_from_residue(res, 'C')
        cb = _atom_xyz_from_residue(res, 'CB')

        if n is not None and ca is not None:
            if i == 0 and c is not None:
                if _atom_xyz_from_residue(res, 'H') is None:
                    _append_atom(additions, 'H', aa_name, resid, place_fourth_atom(c, ca, n, 1.01, 109.5, 120.0))
                if _atom_xyz_from_residue(res, 'H2') is None:
                    _append_atom(additions, 'H2', aa_name, resid, place_fourth_atom(c, ca, n, 1.01, 109.5, -120.0))
            elif original_aa != 'P' and _atom_xyz_from_residue(res, 'H') is None:
                prev = residue_lookup[residues[i - 1]]
                prev_c = _atom_xyz_from_residue(prev, 'C')
                if prev_c is not None:
                    h = _place_opposite_bisector(n, [prev_c, ca], 1.01)
                    if h is not None:
                        _append_atom(additions, 'H', aa_name, resid, h)

        if n is not None and ca is not None and c is not None:
            if original_aa == 'G':
                if _atom_xyz_from_residue(res, 'HA2') is None:
                    _append_atom(additions, 'HA2', aa_name, resid, place_fourth_atom(n, c, ca, 1.09, 109.5, 120.0))
                if _atom_xyz_from_residue(res, 'HA3') is None:
                    _append_atom(additions, 'HA3', aa_name, resid, place_fourth_atom(n, c, ca, 1.09, 109.5, -120.0))
            elif _atom_xyz_from_residue(res, 'HA') is None:
                neighbors = [n, c] if cb is None else [n, c, cb]
                ha = _place_opposite_bisector(ca, neighbors, 1.09)
                if ha is not None:
                    _append_atom(additions, 'HA', aa_name, resid, ha)

    if additions:
        df = pd.concat([df, pd.DataFrame(additions)], ignore_index=True)

    df['_order'] = range(len(df))
    df['_added_h'] = (df['anum'] == 0).astype(int)
    df = df.sort_values(['resid', '_added_h', '_order']).drop(columns=['_order', '_added_h'])
    df = df.reset_index(drop=True)
    df['anum'] = range(1, len(df) + 1)
    return df


def _cb_local_frame(residue_df):
    ca = _atom_xyz_from_residue(residue_df, 'CA')
    cb = _atom_xyz_from_residue(residue_df, 'CB')
    if ca is None or cb is None:
        raise ValueError('Need CA and CB atoms to build a side-chain transplant frame.')

    x_axis = unit_vector(cb - ca, 'CA-CB vector')
    y_axis = None
    for anchor_name in ('N', 'C'):
        anchor = _atom_xyz_from_residue(residue_df, anchor_name)
        if anchor is None:
            continue
        anchor_vec = anchor - ca
        anchor_vec = anchor_vec - np.dot(anchor_vec, x_axis) * x_axis
        if np.linalg.norm(anchor_vec) >= 1e-8:
            y_axis = unit_vector(anchor_vec, f'{anchor_name}-CA plane vector')
            break

    if y_axis is None:
        raise ValueError('Need N or C atom to orient the side-chain transplant frame.')

    z_axis = unit_vector(np.cross(x_axis, y_axis), 'side-chain frame normal')
    y_axis = unit_vector(np.cross(z_axis, x_axis), 'side-chain frame y-axis')
    return np.column_stack([x_axis, y_axis, z_axis])


def Transplant(backbone_df, sequence, rotamer_dir='RotamerLibrary', residue_map=None):
    # Keep the builder backbone/CB atoms, then add side-chain atoms from RotamerLibrary.
    sequence = sequence.upper().strip()
    residue_map = dict(AA1_TO_ROTAMER if residue_map is None else residue_map)
    df = backbone_df.copy()
    sidechain_parts = []

    for resid, aa in enumerate(sequence, start=1):
        residue_name = residue_map[aa]
        df.loc[df['resid'] == resid, 'aa_name'] = residue_name
        keep_atoms = BACKBONE_TEMPLATE_ATOMS | {'CB'}
        df = df[(df['resid'] != resid) | df['atom_name'].isin(keep_atoms)].copy()

        if aa == 'G':
            continue

        target_res = df[df['resid'] == resid].copy()
        template = read_rotamer_template(residue_name, rotamer_dir=rotamer_dir).copy()

        template_ca = _atom_xyz_from_residue(template, 'CA')
        template_cb = _atom_xyz_from_residue(template, 'CB')
        target_ca = _atom_xyz_from_residue(target_res, 'CA')
        target_cb = _atom_xyz_from_residue(target_res, 'CB')
        if template_ca is None or template_cb is None or target_ca is None or target_cb is None:
            raise ValueError(f'Need CA and CB atoms to transplant residue {residue_name} at resid {resid}.')

        template_cb_length = np.linalg.norm(template_cb - template_ca)
        target_cb_direction = unit_vector(target_cb - target_ca, 'builder CA-CB vector')
        placed_cb = target_ca + target_cb_direction * template_cb_length
        cb_mask = (df['resid'] == resid) & (df['atom_name'] == 'CB')
        df.loc[cb_mask, ['x', 'y', 'z']] = placed_cb
        target_res.loc[target_res['atom_name'] == 'CB', ['x', 'y', 'z']] = placed_cb

        rotation = _cb_local_frame(target_res) @ _cb_local_frame(template).T

        side_atoms = template[~template['atom_name'].isin(BACKBONE_TEMPLATE_ATOMS | {'CB'})].copy()
        if side_atoms.empty:
            continue

        coords = side_atoms[['x', 'y', 'z']].to_numpy(dtype=np.float64)
        side_atoms[['x', 'y', 'z']] = (coords - template_cb) @ rotation.T + placed_cb
        side_atoms['resid'] = resid
        side_atoms['aa_name'] = residue_name
        sidechain_parts.append(side_atoms)

    if sidechain_parts:
        df = pd.concat([df] + sidechain_parts, ignore_index=True)

    df['_order'] = range(len(df))
    df = df.sort_values(['resid', '_order']).drop(columns=['_order']).reset_index(drop=True)
    df['anum'] = range(1, len(df) + 1)
    return df


def build_single_peptide_with_local_peptidebuilder(sequence, angles, output_file='temp_sequence.pdb',
                                                   rotamer_dir='RotamerLibrary', residue_map=None):
    # Build one peptide strand from the requested sequence, keeping builder backbone/CB geometry.
    sequence = sequence.upper().strip()
    
    
    scaffold_sequence = sequence
    
    backbone_df = build_target_beta_backbone(scaffold_sequence, angles)
    
    backbone_df = Add_backbone_hydrogen(backbone_df, sequence)
    peptide_df = Transplant(backbone_df, sequence, rotamer_dir=rotamer_dir, residue_map=residue_map)
    _write_simple_pdb(peptide_df, output_file)
    return peptide_df


# In[65]:


def alignment(length, chains, pdbfile, o_filename, order=1):
    pdb=[]
    columns = ['atom', 'anum', 'atom_name', 'aa_name', 'resid', 'x', 'y', 'z']

    with open(pdbfile, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("ATOM"):
            components = line.split()
    
            # Create a dictionary for the current line  
            atom = line[0:6].strip()           # ATOM
            anum = int(line[6:11])            # Atom number
            atom_name = line[11:17].strip()   # Atom name
            aa_name = line[17:21].strip()     # Amino acid name
            resid = int(line[22:26])          # Residue ID
            x = float(line[26:38])            # X-coordinate
            y = float(line[38:46])            # Y-coordinate
            z = float(line[46:54])            # Z-coordinate
    
            row = {
                'atom': atom,
                'anum': anum,
                'atom_name': atom_name,
                'aa_name': aa_name,
                'resid': resid,
                'x': x,
                'y': y,
                'z': z
            }
    
            # Append the dictionary to the list
            pdb.append(row)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(pdb, columns=columns)
    

    if (df[df['resid'] == 1]['aa_name'].iloc[0] != "ACE"):
        filtered_df = df[(df['resid'] <= 10) & (df['atom_name'].isin(['C', 'CA', 'CB']))]
    else:
        filtered_df = df[(df['resid'] > 1) &
                                  (df['resid'] <= 11) &
                                  (df['aa_name'] != "NME") &
                                  (df['atom_name'].isin(['C', 'CA', 'CB']))]

        
    # filtered_df = df[(df['resid'] <= length * chains) & (df['atom_name'].isin(['O', 'C']))]
    data = filtered_df[['x', 'y', 'z']].to_numpy()

    # regular grid covering the domain of the data
    X_min = filtered_df["x"].min()
    X_max = filtered_df["x"].max()
    Y_min = filtered_df["y"].min()
    Y_max = filtered_df["y"].max()
    X,Y = np.meshgrid(np.arange(X_min, X_max, 0.6), np.arange(Y_min, Y_max, 0.5))
    XX = X.flatten()
    YY = Y.flatten()

    ################## fitting sheet plane and get the plane normal vector #################
    # order = 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        # data[:,0] = x-column, data[:,1] = y-column, data[:,2] = z-column
        # 先准备matrix A，column 1 2 分别是 实际的x y数据，column 3只包含 1.
        # 拟合公式 Z=C0*X + C1*Y + C2. Matrix A对应 [X, Y, 1].
        # lstsq 里面,第一个输入matrix A，第二个输入实际Z的值，然后返回拟合公式的坐标值
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]

        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
        print(A)
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)

    # plot points and fitted surface
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
#     ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
#     plt.xlabel('X')
#     plt.ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.axis('equal')
#     ax.axis('tight')
#     plt.show()

    ############calculate the angles between z-axis and the plane normal vector ###################
    nz=[0,0,1]
    # sheets plane normal vector = (C[0], C[1], -1)
    # rotation on x-y plane (rotate by z-axis). v1 = (0,0), *** no need to rotate by z-axis ***
    # v1=[C[0],C[1]]
    # v2=[nz[0],nz[1]]
    # tz = rotation_angle(v1,v2) # rotation from v1 to v2

    # rotation on x-z plane (rotate by y-axis)
    v1=[C[0],-1]
    v2=[nz[0],nz[2]]
    ty = rotation_angle(v1,v2) # rotation from v1 to v2

    # rotation on y-z plane (rotate by x-axis)
    v1=[C[1],-1]
    v2=[nz[1],nz[2]]
    tx = rotation_angle(v1,v2) # rotation from v1 to v2

    ############calculate rotation matrix and rotate the plane normal to z-axis ###################
    Rx_M = np.array([[1, 0, 0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]])
    Ry_M = np.array([[np.cos(ty), 0, np.sin(ty)], [0, 1, 0], [-np.sin(ty), 0, np.cos(ty)]])
    # Rz_M = np.array([[np.cos(ty), 0, np.sin(ty)], [0, 1, 0], [-np.sin(ty), 0, np.cos(ty)]])

    data = df[['x', 'y', 'z']].to_numpy()
    df_temp = df.copy()

    for i, xyz in enumerate(data):
        xyz = Rx_M @ xyz
        xyz = Ry_M @ xyz

        df_temp.at[i, 'x'] = xyz[0]
        df_temp.at[i, 'y'] = xyz[1]
        df_temp.at[i, 'z'] = xyz[2]
    
    ############################# fit each peptide parallel to y-axis and move the sheets to center ###############
    vectors=[]

        # filtered_df = df_temp[(df_temp['resid']>i*length)&(df_temp['resid']<=(i+1)*length)&(df_temp['atom_name'].isin(['C','CA','N']))]
    filtered_df = df_temp[df_temp['atom_name'].isin(['CA'])]
    points = filtered_df[['x', 'y', 'z']].to_numpy()

#       k1, b1, k2, b2 = linear_fitting_3D_points(points)

#         p = 1
#         m = k1 * p
#         n = k2 * p
    m,n,p = linear_fitting(points)
    vector = [m, n, p] # sheet direction vector
    vectors.append(vector)

    vectors=np.array(vectors)
    vector = np.mean(vectors, axis=0)
    
     ############################# align peptide to y-axis and move center to origin ###############
    ny = [0, 1, 0]
    v1=[vector[0],vector[1]]
    v2=[ny[0],ny[1]]
    tz = rotation_angle(v1,v2) # rotation from v1 to v2

    Rz_M = np.array([[np.cos(tz), -np.sin(tz), 0], [np.sin(tz), np.cos(tz), 0], [0, 0, 1]])
    
    data = df_temp[['x', 'y', 'z']].to_numpy()

    for i, xyz in enumerate(data):
        xyz = Rz_M @ xyz
        df_temp.at[i, 'x'] = xyz[0]
        df_temp.at[i, 'y'] = xyz[1]
        df_temp.at[i, 'z'] = xyz[2]
    
    x_sum=0; y_sum=0; z_sum=0
    data_for_center = df_temp[(df_temp['atom_name'].isin(['C', 'CA', 'N'])) &
                              (df_temp['aa_name'] != "ACE") &
                              (df_temp['aa_name'] != "NME")
                             ][['x', 'y', 'z']].to_numpy()
    for i, xyz in enumerate(data_for_center):
        x_sum+=xyz[0]
        y_sum+=xyz[1]
        z_sum+=xyz[2]

    x_sum/=(i+1)
    y_sum/=(i+1)
    z_sum/=(i+1)
    df_temp['x'] = df_temp['x'] - x_sum
    df_temp['y'] = df_temp['y'] - y_sum
    df_temp['z'] = df_temp['z'] - z_sum
    
    ############# rotate even residues toward z-axis and 较大序号基团朝向 y-axis ########################
    if (df[df['resid'] == 1]['aa_name'].iloc[0] != "ACE"):
        res_1_CA_x = df_temp[(df_temp['resid'] == 1) & (df_temp['atom_name'] == 'CA')]['x'].iloc[0]
        res_2_CA_x = df_temp[(df_temp['resid'] == 2) & (df_temp['atom_name'] == 'CA')]['x'].iloc[0]
        res_1_CA_y = df_temp[(df_temp['resid'] == 1) & (df_temp['atom_name'] == 'CA')]['y'].iloc[0]
        res_2_CA_y = df_temp[(df_temp['resid'] == 2) & (df_temp['atom_name'] == 'CA')]['y'].iloc[0]
    else:
        res_1_CA_x = df_temp[(df_temp['resid'] == 2) & (df_temp['atom_name'] == 'CA')]['x'].iloc[0]
        res_2_CA_x = df_temp[(df_temp['resid'] == 3) & (df_temp['atom_name'] == 'CA')]['x'].iloc[0]
        res_1_CA_y = df_temp[(df_temp['resid'] == 2) & (df_temp['atom_name'] == 'CA')]['y'].iloc[0]
        res_2_CA_y = df_temp[(df_temp['resid'] == 3) & (df_temp['atom_name'] == 'CA')]['y'].iloc[0]        

    if(res_1_CA_x > res_2_CA_x):
        print("res_1_CA_x > res_2_CA_x")
        ty = np.pi/2 # 绕 y-axis 逆时针旋转 90°,即绕 y-axis 顺时针旋转90°
    else:
        print("res_1_CA_x < res_2_CA_x")
        ty = -np.pi/2

    df_temp = rotation_coordinates(df_temp, 'y', ty) # 绕 y-axis 旋转 90° or -90

    
    # 判断第一个基团的CA 的y坐标是否大于第二个基团的y坐标,即多肽链C端是否朝向y正方向
    if(res_1_CA_y > res_2_CA_y):
        print("res_1_CA_y > res_2_CA_y")
        tz = np.pi                                        # set 180°
        df_temp = rotation_coordinates(df_temp, 'z', tz)  # 绕 z-axis 旋转 180°


#     df_temp = df_temp[(df_temp['resid'] <= length)]
    
    ############# write out the transformed structure ###################     
    with open(o_filename, 'w') as f:
        for index, row in df_temp.iterrows():
            if len(row['aa_name'])==4:
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])
            else:
                output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                    row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                    row['x'],row['y'],row['z'])

            f.write(output)
            
    return df_temp


# In[66]:


def packpep(single_peptite, chain_length, classes, shifts, chains_per_sheet, format_flag, o_filename, 
            core_residues="e", strand_dist = 4.8, sheet_dist = 11.5, sheet_shift = 2.4):
    """
    single_peptite:   DataFrame, singple peptide PDB file information
    chain_length:     int, number of amino acid of the peptide
    classes:          int, type of steric zipper
    core_residues:    "e" = even number of residues in fibril core (default). "o" = odd number of residues in fibril core
    shifts:           int, move sheet = 1 or = -1 residue relative to another, or 0 = not move
    chains_per_sheet: int, how many peptides per sheet
    """
    # sheet_shift = 2.4 # shift along fibril direction (x in this code)
    # strand_dist = 4.8
    # sheet_dist = 11.5
    
    if (classes == 1):
        class_axis = "x"
        if (core_residues == "e"):
        # sheet cross sectional view

        ###### shifts = 0 ########        ###### shifts = 1 ########        ###### shifts = -1 ########

        # z      7 5 3 1                  # z       7 5 3 1               # z     7 5 3 1 
        # |    C--------N  (SHEET 2)      # |     C--------N  (SHEET 2)   # |   C--------N  (SHEET 2)
        # |     8 6 4 2                   # |      8 6 4 2                # |    8 6 4 2
        # |     2 4 6 8                   # |     2 4 6 8                 # |     2 4 6 8
        # |   N--------C   (SHEET 1)      # |   N--------C    (SHEET 1)   # |   N--------C  (SHEET 1)
        # |    1 3 5 7                    # |    1 3 5 7                  # |    1 3 5 7
        # ----------------> y             # ----------------> y           # ----------------> y             
  
            if (chain_length%2 == 0):
                x_adjust = 0 + sheet_shift
                y_adjust = 2.75 + shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(0/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(5/360)

            elif (chain_length%2 == 1):
                x_adjust = 0 + sheet_shift
                y_adjust = -0.25 + shifts * 3.275
                z_adjust = 0 
                tx_adjust = np.pi*(0/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(5/360)
                
        elif (core_residues == "o"):      
            sheet_dist = -sheet_dist       
        # sheet cross sectional view
        # shift = +1: upper sheets move right, shift = -1: upper sheets move left
        ###### shifts = 0 ########       ###### shifts = 1 ########         ###### shifts = -1 ########

        # z      2 4 6 8                    # z     2 4 6 8                # z     2 4 6 8
        # |    N--------C   (SHEET 1)       # |   N--------C  (SHEET 1)    # |   N--------C    (SHEET 1)
        # |     1 3 5 7                     # |    1 3 5 7                 # |    1 3 5 7
        # |     7 5 3 1                     # |   7 5 3 1                  # |     7 5 3 1
        # |   C--------N    (SHEET 2)       # | C--------N    (SHEET 2)    # |   C--------N    (SHEET 2)
        # |    8 6 4 2                      # |  8 6 4 2                   # |    8 6 4 2
        # ----------------> y               # ----------------> y          # ----------------> y 
            if (chain_length%2 == 0):
                x_adjust = 0 + sheet_shift
                y_adjust = -3.7 - shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(0/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(5/360)

            elif (chain_length%2 == 1):
                x_adjust = 0 + sheet_shift
                y_adjust = -0.25 - shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(0/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(5/360)
                
    elif (classes == 3):
        class_axis = "y"
        if (core_residues == "e"):
        # sheet cross sectional view

        ###### shifts = 0 ########       ###### shifts = 1 ########        ###### shifts = -1 ########

        # z    1 3 5 7                    # z      1 3 5 7                # z   1 3 5 7 
        # |   N--------C  (SHEET 2)       # |    N--------C  (SHEET 2)    # |  N--------C    (SHEET 2)
        # |     2 4 6 8                   # |      2 6 4 8                # |    2 4 6 8
        # |     2 4 6 8                   # |     2 4 6 8                 # |     2 4 6 8
        # |   N--------C  (SHEET 1)       # |   N--------C   (SHEET 1)    # |   N--------C   (SHEET 1)
        # |    1 3 5 7                    # |    1 3 5 7                  # |    1 3 5 7
        # ----------------> y             # ----------------> y           # ----------------> y  
            if (chain_length%2 == 0):
                x_adjust = 0 + sheet_shift
                y_adjust = -0.49 + shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(-2/360)
                ty_adjust = np.pi*(0/360)       
                tz_adjust = np.pi*(4/360)
            elif (chain_length%2 == 1):
                x_adjust = 0 + sheet_shift
                y_adjust = -0.25 + shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(-1/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(4/360)
        # sheet cross sectional view

        ###### shifts = 0 ########       ###### shifts = 1 ########        ###### shifts = -1 ########

        # z     2 4 6 8                   # |      2 4 6 8                 # |     2 4 6 8
        # |   N--------C  (SHEET 1)       # |    N--------C  (SHEET 1)     # |   N--------C   (SHEET 1)
        # |    1 3 5 7                    # |     1 3 5 7                  # |    1 3 5 7              
        # |    1 3 5 7                    # |    1 3 5 7                   # |     1 3 5 7
        # |   N--------C  (SHEET 2)       # |  N--------C    (SHEET 2)     # |    N--------C  (SHEET 2)
        # |     2 4 6 8                   # |   2 6 4 8                    # |      2 6 4 8
        # ----------------> y             # ----------------> y            # ----------------> y              
        elif (core_residues == "o"):
            sheet_dist = -sheet_dist
            if (chain_length%2 == 0):
                x_adjust = 0 + sheet_shift
                y_adjust = -0.49 - shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(-2/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(4/360)
            elif (chain_length%2 == 1):
                x_adjust = 0 + sheet_shift
                y_adjust = -0.25 - shifts * 3.275
                z_adjust = 0
                tx_adjust = np.pi*(-1/360)
                ty_adjust = np.pi*(0/360)
                tz_adjust = np.pi*(4/360)
                
    elif (classes == 4):
    # sheet cross sectional view

    ###### shifts = 0 ########       ###### shifts = 1 ########        ###### shifts = -1 ########

    # z    8 6 4 2                    # |      8 6 4 2                 # |   8 6 4 2
    # |   C--------N   (SHEET 2)      # |    C--------N  (SHEET 2)     # |  C--------N    (SHEET 2)
    # |     7 5 3 1                   # |     7 5 3 1                  # |    7 5 3 1              
    # |     2 4 6 8                   # |    2 4 6 8                   # |     2 4 6 8
    # |   N--------C   (SHEET 1)      # |  N--------C    (SHEET 1)     # |   N--------C   (SHEET 1)
    # |    1 3 5 7                    # |   1 3 5 7                    # |    1 3 5 7              
    # ----------------> y             # ----------------> y            # ----------------> y          
        class_axis = "z"
        if (chain_length%2 == 0):
            x_adjust = 0 + sheet_shift
            y_adjust = -0.48 + shifts * 3.275
            z_adjust = 0
            tx_adjust = np.pi*(-1/180)
            ty_adjust = np.pi*(0/180)            
            tz_adjust = np.pi*(0/180)
        elif (chain_length%2 == 1):
            x_adjust = 0 + sheet_shift
            y_adjust = -0.25 + shifts * 3.275
            z_adjust = 0
            tx_adjust = np.pi*(-0.5/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)
    
    elif (classes == 2):
    # sheet cross sectional view

    ###### shifts = 0 ########       ###### shifts = 1 ########        ###### shifts = -1 ########

    # z      2 4 6 8                  # |      2 4 6 8                 # |     2 4 6 8
    # |    N--------C  (SHEET 2)      # |    N--------C  (SHEET 2)     # |   N--------C   (SHEET 2)
    # |     1 3 5 7                   # |     1 3 5 7                  # |    1 3 5 7              
    # |     2 4 6 8                   # |    2 4 6 8                   # |     2 4 6 8
    # |   N--------C   (SHEET 1)      # |  N--------C    (SHEET 1)     # |   N--------C   (SHEET 1)
    # |    1 3 5 7                    # |   1 3 5 7                    # |    1 3 5 7              
    # ----------------> y             # ----------------> y            # ----------------> y          
        class_axis = "none"

        if (chain_length%2 == 0):
            x_adjust = 0 + sheet_shift
            y_adjust = 2.75 + shifts * 3.275
            z_adjust = 0
            tx_adjust = np.pi*(-1/180)
            ty_adjust = np.pi*(0/180)            
            tz_adjust = np.pi*(0/180)
        elif (chain_length%2 == 1):
            x_adjust = 0 + sheet_shift
            y_adjust = 2.75 + shifts * 3.275
            z_adjust = 0
            tx_adjust = np.pi*(-1/180)
            ty_adjust = np.pi*(0/180)            
            tz_adjust = np.pi*(0/180)

###### set parameter ###############

    sheets=2
    anum_max = single_peptite['anum'].max()
    resid_max = single_peptite['resid'].max()
    peptides = pd.DataFrame([]) # empty dataframe
    # print(anum_max)
###### pack first sheet ##################
    for chain in range(chains_per_sheet):
        new_peptide = single_peptite.copy()
        new_peptide['x']     = new_peptide['x']     + chain * strand_dist
        new_peptide['anum']  = new_peptide['anum']  + chain * anum_max
        new_peptide['resid'] = new_peptide['resid'] + chain * resid_max
        peptides = pd.concat([peptides, new_peptide])
        
###### create a temp peptide, rotate to form class-1 sheets and move upward by sheets-ditance, the even residues are facing each other ####
    temp_pep = single_peptite.copy()                            # copy a original peptide
    temp_pep['anum']  = temp_pep['anum'] + chains_per_sheet * anum_max             # shift atom number 
    temp_pep['resid'] = temp_pep['resid'] + chains_per_sheet * resid_max   # shift residue number
    
    
    theta = -np.pi                                              # construct 2nd peptide model to fit different classes
    temp_pep = rotation_coordinates(temp_pep, class_axis, theta)
    
    temp_pep['z'] = temp_pep['z'] + sheet_dist                  # fine adjustment
    
    temp_pep['x'] = temp_pep['x'] + x_adjust
    temp_pep['y'] = temp_pep['y'] + y_adjust
    temp_pep['z'] = temp_pep['z'] + z_adjust
    temp_pep = rotation_coordinates(temp_pep, 'x', tx_adjust)
    temp_pep = rotation_coordinates(temp_pep, 'y', ty_adjust)
    temp_pep = rotation_coordinates(temp_pep, 'z', tz_adjust)

    
###### pack second sheet ####################
    for chain in range(chains_per_sheet):
        new_peptide = temp_pep.copy()
        new_peptide['x'] = new_peptide['x'] + chain * strand_dist
        new_peptide['anum'] = new_peptide['anum'] + chain * anum_max
        new_peptide['resid'] = new_peptide['resid'] + chain * resid_max
        peptides = pd.concat([peptides, new_peptide])



    if format_flag == 0:
        last_resid = 1
        with open(o_filename, 'w') as f:
            for index, row in peptides.iterrows():
                resid = int(row['resid'])
                if len(row['aa_name'])==4:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                        row['x'],row['y'],row['z'])
                else:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                        row['x'],row['y'],row['z'])
    
                f.write(output)
                last_resid = resid
    
    elif format_flag == 1:
        last_resid = 1
        with open(o_filename, 'w') as f:
            for index, row in peptides.iterrows():
                resid = int(row['resid'])
                
                if len(row['aa_name'])==4:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'][1:],row['resid'],
                        row['x'],row['y'],row['z'])
                else:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                        row['x'],row['y'],row['z'])
    
                if (last_resid % chain_length == 0) & (resid % chain_length ==1):
                     f.write('TER\n')
                        
                f.write(output)
                last_resid = resid
            
def packpep_antiparallel(single_peptite, chain_length, classes, shifts, chains_per_sheet, 
                         format_flag, o_filename, core_residues="e", strand_dist = 4.8, sheet_dist = 11.5, sheet_shift = 0):
    # strand_dist = 4.8
    # sheet_dist = 11.5
    # sheet_shift = 2.4
    
    if (classes == 5):
        class_axis = "x"
        if (core_residues == "e"):
        # sheet cross sectional view

        ###### shifts = 0 ########                   ###### shifts = 1 ########                         
        # Z                                          # z                                     
        # |              1 3 5 7                     # |                1 3 5 7                      
        # |     7 5 3 1 N--------C                   # |       7 5 3 1 N--------C                    
        # |   C--------N  2 4 6 8     (SHEET 2)      # |     C--------N  2 4 6 8   (SHEET 2)         
        # |    8 6 4 2    8 6 4 2                    # |      8 6 4 2   8 6 4 2                      
        # |    2 4 6 8   C--------N                  # |     2 4 6 8   C--------N                    
        # |  N--------C    7 5 3 1    (SHEET 1)      # |   N--------C    7 5 3 1   (SHEET 1)         
        # |   1 3 5 7                                # |    1 3 5 7                  
        # ---------------------------> y             # ----------------------------> y                        
  
        ###### shifts = -1 ########                            
        # Z                                                        
        # |            1 3 5 7                                   
        # |   7 5 3 1 N--------C                          
        # | C--------N  2 4 6 8       (SHEET 2)               
        # |  8 6 4 2   8 6 4 2                                      
        # |   2 4 6 8 C--------N                                 
        # | N--------C  7 5 3 1       (SHEET 1)            
        # |  1 3 5 7                                              
        # ---------------------------> y 
    
            if (chain_length%2 == 0):
                x_adjust = 0 + sheet_shift
                y_adjust = 3.3 + shifts * 3.461
                z_adjust = 0
                tx_adjust = np.pi*(0/180)
                ty_adjust = np.pi*(0/180)
                tz_adjust = np.pi*(0/180)
            
            elif (chain_length%2 == 1):
                x_adjust = 0 + sheet_shift
                y_adjust = 0 + shifts * 3.461
                z_adjust = 0
                tx_adjust = np.pi*(0/180)
                ty_adjust = np.pi*(0/180)
                tz_adjust = np.pi*(0/180)                   
                
        if (core_residues == "o"):
            sheet_dist = -sheet_dist
        # sheet cross sectional view

        ###### shifts = 0 ########                   ###### shifts = 1 ########                         
        # Z                                          # z                                      
        # |               2 4 6 8                    # |                 2 4 6 8                      
        # |    8 6 4 2  N--------C                   # |      8 6 4 2  N--------C                    
        # |   C--------N 1 3 5 7      (SHEET 1)      # |     C--------N 1 3 5 7     (SHEET 1)         
        # |     7 5 3 1   7 5 3 1                    # |       7 5 3 1    7 5 3 1                     
        # |   1 3 5 7   C--------N                   # |    1 3 5 7    C--------N                    
        # |  N--------C  8 6 4 2      (SHEET 2)      # |   N--------C   8 6 4 2     (SHEET 2)         
        # |    2 4 6 8                               # |     2 4 6 8                  
        # ---------------------------> y             # ----------------------------> y                        
  
        ###### shifts = -1 ########                                           
        # Z                                                                               
        # |              2 4 6 8                                        
        # |   8 6 4 2  N--------C                                     
        # |  C--------N 1 3 5 7      (SHEET 1)            
        # |    7 5 3 1   7 5 3 1                                   
        # |   1 3 5 7  C--------N                         
        # |  N--------C 8 6 4 2      (SHEET 2)        
        # |    2 4 6 8                                               
        # ---------------------------> y  
    
            if (chain_length%2 == 0):
                x_adjust = 0 + sheet_shift
                y_adjust = 3.3 + shifts * 3.461
                z_adjust = 0
                tx_adjust = np.pi*(0/180)
                ty_adjust = np.pi*(0/180)
                tz_adjust = np.pi*(0/180)

            elif (chain_length%2 == 1):
                x_adjust = 0 + sheet_shift
                y_adjust = shifts * 3.461
                z_adjust = 0
                tx_adjust = np.pi*(0/180)
                ty_adjust = np.pi*(0/180)
                tz_adjust = np.pi*(0/180)
        
    ###### set parameter ###############
        sheets=2
        anum_max = single_peptite['anum'].max()
        resid_max = single_peptite['resid'].max()

        anti_pep = single_peptite.copy()                     # construct antiparallel peptide in same sheet
        theta = np.pi*(180/360)                                       
        anti_pep = rotation_coordinates(anti_pep, "z", theta)

    ############################# align anti-peptide to single-pep ###############
        filtered_df = anti_pep[(anti_pep['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        anti_pep_vector = np.array([m, n, p])
        if (anti_pep_vector[1] > 0):  # If y-component is positive, flip direction
            anti_pep_vector = -anti_pep_vector

        filtered_df = single_peptite[(single_peptite['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        singple_pep_vector = np.array([m, n, p])
        if (singple_pep_vector[1] > 0):  # If y-component is positive, flip direction
            singple_pep_vector = -singple_pep_vector

        v1=[anti_pep_vector[0],anti_pep_vector[1]]
        v2=[singple_pep_vector[0],singple_pep_vector[1]]
        tz = rotation_angle(v1,v2) # rotation from v1 to v2
        anti_pep = rotation_coordinates(anti_pep, 'z', tz)
        
        first_resid = single_peptite[single_peptite['atom_name'] == 'CA']['resid'].min()
        ca1_y = single_peptite[(single_peptite['atom_name'] == 'CA') & (single_peptite['resid'] == first_resid)].iloc[0]['y']
        last_resid = anti_pep[anti_pep['atom_name'] == 'CA']['resid'].max()
        ca2_y = anti_pep[(anti_pep['atom_name'] == 'CA') & (anti_pep['resid'] == last_resid)].iloc[0]['y']

        
        
        anti_pep_adjust = ca1_y - ca2_y 
        print(first_resid, last_resid, "ca1_y", ca1_y, "ca2_y",ca2_y, anti_pep_adjust)    
        if (chain_length%2 == 0):
            anti_pep_adjust += 3.461
    ################### anti_pep adjustment ###############  
        anti_pep['y']   = anti_pep['y'] + anti_pep_adjust
        
    ###### pack first sheet ##################
        peptides = pd.DataFrame([]) # empty dataframe
        for chain in range(chains_per_sheet):
            if (chain % 2 == 0):
                new_peptide = single_peptite.copy()
            elif (chain % 2 == 1):
                new_peptide = anti_pep.copy()

            new_peptide['x']     = new_peptide['x']     + chain * strand_dist
            new_peptide['anum']  = new_peptide['anum']  + chain * anum_max
            new_peptide['resid'] = new_peptide['resid'] + chain * resid_max
            peptides = pd.concat([peptides, new_peptide])
            
    ###### build peptide for second sheet
        temp_sheet = peptides.copy()                                            # copy a original peptide
        temp_sheet['anum']  = temp_sheet['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_sheet['resid'] = temp_sheet['resid'] + chains_per_sheet * resid_max   # shift residue number
        theta = -np.pi                                                         # build pep for 2nd sheet
        temp_sheet = rotation_coordinates(temp_sheet, "x", theta)

        temp_sheet['z'] = temp_sheet['z'] + sheet_dist
        temp_sheet['x'] = temp_sheet['x'] + x_adjust
        temp_sheet['y'] = temp_sheet['y'] + y_adjust
        temp_sheet['z'] = temp_sheet['z'] + z_adjust
        temp_sheet = rotation_coordinates(temp_sheet,'x',tx_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'y',ty_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'z',tz_adjust)
        
    ###### pack second sheet ####################
        peptides = pd.concat([peptides, temp_sheet])

        
    elif (classes == 6):
        class_axis = "none"

        # sheet cross sectional view

        ###### shifts = 0 ########                   ###### shifts = 1 ########                         
        # Z                                          # z                                     
        # |              8 6 4 2                     # |                8 6 4 2                      
        # |     2 4 6 8 C--------N                   # |       2 4 6 8 C--------N                    
        # |   N--------C  7 5 3 1     (SHEET 2)      # |     N--------C  7 5 3 1   (SHEET 2)         
        # |    1 3 5 7    8 6 4 2                    # |      1 3 5 7   8 6 4 2                      
        # |    2 4 6 8   C--------N                  # |     2 4 6 8   C--------N                    
        # |  N--------C    7 5 3 1    (SHEET 1)      # |   N--------C    7 5 3 1   (SHEET 1)         
        # |   1 3 5 7                                # |    1 3 5 7                  
        # ---------------------------> y             # ----------------------------> y                        
  
        ###### shifts = -1 ########                            
        # Z                                                        
        # |            8 6 4 2                                   
        # |   2 4 6 8 C--------N                          
        # | N--------C  7 5 3 1       (SHEET 2)               
        # |  1 3 5 7   8 6 4 2                                      
        # |   2 4 6 8 C--------N                                 
        # | N--------C  7 5 3 1       (SHEET 1)            
        # |  1 3 5 7                                              
        # ---------------------------> y 
    
        if (chain_length%2 == 0):
            x_adjust = 0 + sheet_shift
            y_adjust = 3.461 + shifts * 3.461
            z_adjust = 0
            tx_adjust = np.pi*(0/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)
            
        elif (chain_length%2 == 1):
            x_adjust = 0 + sheet_shift
            y_adjust = 3.461 + shifts * 3.461
            z_adjust = 0
            tx_adjust = np.pi*(0/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)                
            
    ###### set parameter ###############
        sheets=2
        anum_max = single_peptite['anum'].max()
        resid_max = single_peptite['resid'].max()

        anti_pep = single_peptite.copy()                     # construct antiparallel peptide in same sheet
        theta = np.pi                               
        anti_pep = rotation_coordinates(anti_pep, "z", theta)

    ############################# align anti-peptide to single-pep ###############
        filtered_df = anti_pep[(anti_pep['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        anti_pep_vector = np.array([m, n, p])
        if (anti_pep_vector[1] > 0):  # If y-component is positive, flip direction
            anti_pep_vector = -anti_pep_vector

        filtered_df = single_peptite[(single_peptite['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        singple_pep_vector = np.array([m, n, p])
        if (singple_pep_vector[1] > 0):  # If y-component is positive, flip direction
            singple_pep_vector = -singple_pep_vector

        v1=[anti_pep_vector[0],anti_pep_vector[1]]
        v2=[singple_pep_vector[0],singple_pep_vector[1]]
        tz = rotation_angle(v1,v2) # rotation from v1 to v2
        anti_pep = rotation_coordinates(anti_pep, 'z', tz)
        
        first_resid = single_peptite[single_peptite['atom_name'] == 'CA']['resid'].min()
        ca1_y = single_peptite[(single_peptite['atom_name'] == 'CA') & (single_peptite['resid'] == first_resid)].iloc[0]['y']
        last_resid = anti_pep[anti_pep['atom_name'] == 'CA']['resid'].max()
        ca2_y = anti_pep[(anti_pep['atom_name'] == 'CA') & (anti_pep['resid'] == last_resid)].iloc[0]['y']
        
        anti_pep_adjust = ca1_y - ca2_y
        if (chain_length%2 == 0):
            anti_pep_adjust += 3.461
    ################### anti_pep adjustment ###############  
        anti_pep['y']   = anti_pep['y'] + anti_pep_adjust

        peptides = pd.DataFrame([]) # empty dataframe

    ###### pack first sheet ##################
        for chain in range(chains_per_sheet):
            if (chain % 2 == 0):
                new_peptide = single_peptite.copy()
            elif (chain % 2 == 1):
                new_peptide = anti_pep.copy()

            new_peptide['x']     = new_peptide['x']     + chain * strand_dist
            new_peptide['anum']  = new_peptide['anum']  + chain * anum_max
            new_peptide['resid'] = new_peptide['resid'] + chain * resid_max
            peptides = pd.concat([peptides, new_peptide])
    ###### build peptide for second sheet        
        temp_pep = single_peptite.copy()                                       # copy a original peptide
        temp_pep['anum']  = temp_pep['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_pep['resid'] = temp_pep['resid'] + chains_per_sheet * resid_max   # shift residue number

        temp_anti_pep = anti_pep.copy()                                        # build anti-pep in 2nd sheet
        temp_anti_pep['anum']  = temp_anti_pep['anum'] * anum_max                        # shift atom number 
        temp_anti_pep['resid'] = temp_anti_pep['resid'] + chains_per_sheet * resid_max   # shift residue number

    ###### build peptide for second sheet
        temp_sheet = peptides.copy()                                            # copy a original peptide
        temp_sheet['anum']  = temp_sheet['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_sheet['resid'] = temp_sheet['resid'] + chains_per_sheet * resid_max   # shift residue number

        temp_sheet['z'] = temp_sheet['z'] + sheet_dist
        temp_sheet['x'] = temp_sheet['x'] + x_adjust
        temp_sheet['y'] = temp_sheet['y'] + y_adjust
        temp_sheet['z'] = temp_sheet['z'] + z_adjust
        temp_sheet = rotation_coordinates(temp_sheet,'x',tx_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'y',ty_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'z',tz_adjust)
        
    ###### pack second sheet ####################
        peptides = pd.concat([peptides, temp_sheet])
    
    elif (classes == 7):
        class_axis = "none"

        # sheet cross sectional view

        ###### shifts = 0 ########                   ###### shifts = 1 ########                         
        # Z                                          # z                                     
        # |               7 5 3 1                    # |                  7 5 3 1                      
        # |     2 4 6 8 C--------N                   # |       2 4 6 8  C--------N                    
        # |   N--------C 8 6 4 2     (SHEET 2)       # |     N--------C  8 6 4 2   (SHEET 2)         
        # |    1 3 5 7   7 5 3 1                     # |      1 3 5 7   7 5 3 1                      
        # |    2 4 6 8 C--------N                    # |     2 4 6 8  C--------N                    
        # |  N--------C 8 6 4 2    (SHEET 1)         # |   N--------C  8 6 4 2     (SHEET 1)         
        # |   1 3 5 7                                # |    1 3 5 7                  
        # ---------------------------> y             # ----------------------------> y                        
  
        ###### shifts = -1 ########                            
        # Z                                                        
        # |             7 5 3 1                                   
        # |   2 4 6 8 C--------N                          
        # | N--------C 8 6 4 2        (SHEET 2)               
        # |  1 3 5 7    7 5 3 1                                      
        # |   2 4 6 8 C--------N                                 
        # | N--------C 8 6 4 2        (SHEET 1)            
        # |  1 3 5 7                                              
        # ---------------------------> y 
    
        if (chain_length%2 == 0):
            x_adjust = 0 + sheet_shift
            y_adjust = 3.461 + shifts * 3.461
            z_adjust = 0
            tx_adjust = np.pi*(0/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)
            
        elif (chain_length%2 == 1):
            x_adjust = 0 + sheet_shift
            y_adjust = 3.461 + shifts * 3.461
            z_adjust = 0
            tx_adjust = np.pi*(0/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)                
            
    ###### set parameter ###############
        sheets=2
        anum_max = single_peptite['anum'].max()
        resid_max = single_peptite['resid'].max()

        anti_pep = single_peptite.copy()                     # construct antiparallel peptide in same sheet
        theta = np.pi                                      
        anti_pep = rotation_coordinates(anti_pep, "x", theta)

    ############################# align anti-peptide to single-pep ###############
        filtered_df = anti_pep[(anti_pep['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        anti_pep_vector = np.array([m, n, p])
        if (anti_pep_vector[1] > 0):  # If y-component is positive, flip direction
            anti_pep_vector = -anti_pep_vector

        filtered_df = single_peptite[(single_peptite['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        singple_pep_vector = np.array([m, n, p])
        if (singple_pep_vector[1] > 0):  # If y-component is positive, flip direction
            singple_pep_vector = -singple_pep_vector

        v1=[anti_pep_vector[0],anti_pep_vector[1]]
        v2=[singple_pep_vector[0],singple_pep_vector[1]]
        tz = rotation_angle(v1,v2) # rotation from v1 to v2
        anti_pep = rotation_coordinates(anti_pep, 'z', tz)
        
#         ca1_y = single_peptite[(single_peptite['atom_name'] == 'CA') & (single_peptite['resid'] == 1)].iloc[0][['y']]
#         ca2_y = anti_pep[(anti_pep['atom_name'] == 'CA') & (anti_pep['resid'] == resid_max)].iloc[0][['y']]
#         anti_pep_adjust = ca2_y.values[0] - ca1_y.values[0] 

        first_resid = single_peptite[single_peptite['atom_name'] == 'CA']['resid'].min()
        ca1_y = single_peptite[(single_peptite['atom_name'] == 'CA') & (single_peptite['resid'] == first_resid)].iloc[0]['y']
        last_resid = anti_pep[anti_pep['atom_name'] == 'CA']['resid'].max()
        ca2_y = anti_pep[(anti_pep['atom_name'] == 'CA') & (anti_pep['resid'] == last_resid)].iloc[0]['y']

        anti_pep_adjust = ca1_y - ca2_y

        if (chain_length%2 == 1):
            anti_pep_adjust += 3.461
    ################### anti_pep adjustment ###############  
        anti_pep['y']   = anti_pep['y'] + anti_pep_adjust

        peptides = pd.DataFrame([]) # empty dataframe

    ###### pack first sheet ##################
        for chain in range(chains_per_sheet):
            if (chain % 2 == 0):
                new_peptide = single_peptite.copy()
            elif (chain % 2 == 1):
                new_peptide = anti_pep.copy()

            new_peptide['x']     = new_peptide['x']     + chain * strand_dist
            new_peptide['anum']  = new_peptide['anum']  + chain * anum_max
            new_peptide['resid'] = new_peptide['resid'] + chain * resid_max
            peptides = pd.concat([peptides, new_peptide])
    ###### build peptide for second sheet        
        temp_pep = single_peptite.copy()                                       # copy a original peptide
        temp_pep['anum']  = temp_pep['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_pep['resid'] = temp_pep['resid'] + chains_per_sheet * resid_max   # shift residue number

        temp_anti_pep = anti_pep.copy()                                        # build anti-pep in 2nd sheet
        temp_anti_pep['anum']  = temp_anti_pep['anum'] * anum_max                        # shift atom number 
        temp_anti_pep['resid'] = temp_anti_pep['resid'] + chains_per_sheet * resid_max   # shift residue number

    ###### build peptide for second sheet
        temp_sheet = peptides.copy()                                            # copy a original peptide
        temp_sheet['anum']  = temp_sheet['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_sheet['resid'] = temp_sheet['resid'] + chains_per_sheet * resid_max   # shift residue number

        temp_sheet['z'] = temp_sheet['z'] + sheet_dist
        temp_sheet['x'] = temp_sheet['x'] + x_adjust
        temp_sheet['y'] = temp_sheet['y'] + y_adjust
        temp_sheet['z'] = temp_sheet['z'] + z_adjust
        temp_sheet = rotation_coordinates(temp_sheet,'x',tx_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'y',ty_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'z',tz_adjust)
        
    ###### pack second sheet ####################
        peptides = pd.concat([peptides, temp_sheet])   
    elif (classes == 8):
        class_axis = "z"

        # sheet cross sectional view

        ###### shifts = 0 ########                   ###### shifts = 1 ########                         
        # Z                                          # z                                     
        # |             8 6 4 2                      # |               8 6 4 2                      
        # |   1 3 5 7  C--------N                    # |     1 3 5 7  C--------N                    
        # |  N--------C  7 5 3 1    (SHEET 2)        # |    N--------C  7 5 3 1   (SHEET 2)         
        # |    2 4 6 8   7 5 3 1                     # |      2 4 6 8  7 5 3 1                      
        # |    2 4 6 8 C--------N                    # |     2 4 6 8 C--------N                    
        # |  N--------C 8 6 4 2      (SHEET 1)       # |   N--------C 8 6 4 2     (SHEET 1)         
        # |   1 3 5 7                                # |    1 3 5 7                  
        # ---------------------------> y             # ----------------------------> y                        
  
        ###### shifts = -1 ########                            
        # Z                                                        
        # |            8 6 4 2                                   
        # |  1 3 5 7  C--------N                          
        # | N--------C  7 5 3 1        (SHEET 2)               
        # |   2 4 6 8    7 5 3 1                                      
        # |    2 4 6 8 C--------N                                 
        # |  N--------C 8 6 4 2        (SHEET 1)            
        # |   1 3 5 7                                              
        # ---------------------------> y 
    
        if (chain_length%2 == 0):
            x_adjust = 0 + sheet_shift
            y_adjust = 0.2 + shifts * 3.461
            z_adjust = 0
            tx_adjust = np.pi*(0/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)
            
        elif (chain_length%2 == 1):
            x_adjust = 0 + sheet_shift
            y_adjust = 3.461 + shifts * 3.461
            z_adjust = 0
            tx_adjust = np.pi*(0/180)
            ty_adjust = np.pi*(0/180)
            tz_adjust = np.pi*(0/180)                
            
    ###### set parameter ###############
        sheets=2
        anum_max = single_peptite['anum'].max()
        resid_max = single_peptite['resid'].max()

        anti_pep = single_peptite.copy()                     # construct antiparallel peptide in same sheet
        theta = np.pi                                      
        anti_pep = rotation_coordinates(anti_pep, "x", theta)

    ############################# align anti-peptide to single-pep ###############
        filtered_df = anti_pep[(anti_pep['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        anti_pep_vector = np.array([m, n, p])
        if (anti_pep_vector[1] > 0):  # If y-component is positive, flip direction
            anti_pep_vector = -anti_pep_vector

        filtered_df = single_peptite[(single_peptite['atom_name'].isin(['CA']))]
        points = filtered_df[['x', 'y', 'z']].to_numpy()
        m,n,p = linear_fitting(points)
        singple_pep_vector = np.array([m, n, p])
        if (singple_pep_vector[1] > 0):  # If y-component is positive, flip direction
            singple_pep_vector = -singple_pep_vector

        v1=[anti_pep_vector[0],anti_pep_vector[1]]
        v2=[singple_pep_vector[0],singple_pep_vector[1]]
        tz = rotation_angle(v1,v2) # rotation from v1 to v2
        anti_pep = rotation_coordinates(anti_pep, 'z', tz)
        
        first_resid = single_peptite[single_peptite['atom_name'] == 'CA']['resid'].min()
        ca1_y = single_peptite[(single_peptite['atom_name'] == 'CA') & (single_peptite['resid'] == first_resid)].iloc[0]['y']
        last_resid = anti_pep[anti_pep['atom_name'] == 'CA']['resid'].max()
        ca2_y = anti_pep[(anti_pep['atom_name'] == 'CA') & (anti_pep['resid'] == last_resid)].iloc[0]['y']
        
        anti_pep_adjust = ca1_y - ca2_y
        if (chain_length%2 == 1):
            anti_pep_adjust += 3.461
    ################### anti_pep adjustment ###############  
        anti_pep['y']   = anti_pep['y'] + anti_pep_adjust

        peptides = pd.DataFrame([]) # empty dataframe

    ###### pack first sheet ##################
        for chain in range(chains_per_sheet):
            if (chain % 2 == 0):
                new_peptide = single_peptite.copy()
            elif (chain % 2 == 1):
                new_peptide = anti_pep.copy()

            new_peptide['x']     = new_peptide['x']     + chain * strand_dist
            new_peptide['anum']  = new_peptide['anum']  + chain * anum_max
            new_peptide['resid'] = new_peptide['resid'] + chain * resid_max
            peptides = pd.concat([peptides, new_peptide])
    ###### build peptide for second sheet        
        temp_pep = single_peptite.copy()                                       # copy a original peptide
        temp_pep['anum']  = temp_pep['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_pep['resid'] = temp_pep['resid'] + chains_per_sheet * resid_max   # shift residue number

        temp_anti_pep = anti_pep.copy()                                        # build anti-pep in 2nd sheet
        temp_anti_pep['anum']  = temp_anti_pep['anum'] * anum_max                        # shift atom number 
        temp_anti_pep['resid'] = temp_anti_pep['resid'] + chains_per_sheet * resid_max   # shift residue number

    ###### build peptide for second sheet
        temp_sheet = peptides.copy()                                            # copy a original peptide
        temp_sheet['anum']  = temp_sheet['anum'] + chains_per_sheet * anum_max                        # shift atom number 
        temp_sheet['resid'] = temp_sheet['resid'] + chains_per_sheet * resid_max   # shift residue number

        theta = -np.pi                                                         # build pep for 2nd sheet
        temp_sheet = rotation_coordinates(temp_sheet, class_axis, theta)
        
        temp_sheet['z'] = temp_sheet['z'] + sheet_dist
        temp_sheet['x'] = temp_sheet['x'] + x_adjust + (chains_per_sheet-1) * strand_dist
        temp_sheet['y'] = temp_sheet['y'] + y_adjust
        temp_sheet['z'] = temp_sheet['z'] + z_adjust
        temp_sheet = rotation_coordinates(temp_sheet,'x',tx_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'y',ty_adjust)
        temp_sheet = rotation_coordinates(temp_sheet,'z',tz_adjust)
        
    ###### pack second sheet ####################
        peptides = pd.concat([peptides, temp_sheet])    
    

    ###### output sheets ####################
    if format_flag == 0:
        last_resid = 1
        with open(o_filename, 'w') as f:
            for index, row in peptides.iterrows():
                if len(row['aa_name'])==4:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                        row['x'],row['y'],row['z'])
                else:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                        row['x'],row['y'],row['z'])
                f.write(output)
    
    elif format_flag == 1:
        last_resid = 1
        with open(o_filename, 'w') as f:
            for index, row in peptides.iterrows():
                resid = int(row['resid'])
                
                if len(row['aa_name'])==4:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'][1:],row['resid'],
                        row['x'],row['y'],row['z'])
                else:
                    output= "ATOM{:7d}{:^6}{:4}{:5d}{:12.3f}{:8.3f}{:8.3f}\n".format(
                        row['anum'],row['atom_name'],row['aa_name'],row['resid'],
                        row['x'],row['y'],row['z'])
    
                if (last_resid % chain_length == 0) & (resid % chain_length ==1):
                     f.write('TER\n')
                        
                f.write(output)
                last_resid = resid       

 


# In[67]:


def main():
    parser = argparse.ArgumentParser(description="Peptide builder")

    parser.add_argument('-seq', type=str, required=True, help="Peptide sequence (e.g., GNNQQNY)")
    parser.add_argument('-c', type=int, required=True, choices=range(1, 9), help="Class (1-8 only)")
    parser.add_argument('-sh', type=float, default=0.0, help="Shift value: -1, 0, or 1")
    parser.add_argument('-n', type=int, default=8, help="Number of chains (default: 8)")
    parser.add_argument('-p', type=int, default=0, help="terminal patch uncap = 0 (default), NME capped=1, NHE cappped=2")
    parser.add_argument('-r', type=str, default="e", help="residue number to be packed in the sheets even = e (default), odd = o")
    parser.add_argument('-f', type=int, default=0, help="format flag Pep AD = 0 (default), AMBER = 1")
    parser.add_argument('-d1', type=float, default=4.8, help="strand-strand distance (default=4.8)")
    parser.add_argument('-d2', type=float, default=11.5, help="sheet-sheet distance (default=11.5)")
    parser.add_argument('-d3', type=float, default=0, help="sheet2 moves along fibril axis (default = 0)")
    parser.add_argument('-o', type=str, default="default", help="output file name")
    
    args = parser.parse_args()

    print("Sequence:", args.seq)
    print("Class:", args.c)
    print("Shift:", args.sh)
    print("Number of chains:", args.n)
    print("Cap flag:", args.p)
   
    sequence = args.seq
    classes  = args.c
    shift    = args.sh
    chains   = args.n 
    cap_flag = args.p
    residue_num = args.r
    format_flag = args.f
    strand_dist = args.d1
    sheet_dist = args.d2
    sheet_shift = args.d3
    final_filename  = args.o
    
    if final_filename == "default":
        final_filename = "Class{}_{}mer.pdb".format(classes,residue_num)
    elif final_filename[-4:].lower() == ".pdb":
        pass
    else:
        final_filename = final_filename + ".pdb"
    
    if classes in [1, 2, 3, 4]: #parallel sheets
        angles = [-119, 113]
    elif classes in [5, 6, 7, 8]:
        angles = [-139, 135]
    else:
        raise SystemExit("Cannot determine Class. Stop.")

    ###########################################
    ########## create single pep ##############
    ###########################################
    build_single_peptide_with_local_peptidebuilder(sequence, angles, output_file="temp_pep.pdb")

    ###########################################
    ############# Modify termini ##############
    ###########################################
    if cap_flag == 0:      # Add hydrogen at N-terminus and remove H at C-terminus
        add_NH_remove_OH("temp_pep.pdb", "pep.pdb")
        length=len(sequence)
        
    elif cap_flag in [1, 2]:    # Add Caps to peptide
        add_caps("temp_pep.pdb", "pep.pdb", cap_flag)
        length=len(sequence) + 2
    else:
        raise SystemExit("Cannot determine caps or not. Stop")
        
    ###########################################
    ############# PEP ALIGNMENT ###############
    ###########################################
    single_pep = alignment(length,1,"pep.pdb", "pep_aligned.pdb", order = 1)
    
    ###########################################
    ############## Pack sheets ################
    ###########################################
    if classes in [1, 2, 3, 4]: #parallel sheets
        packpep(single_pep, length, classes, shift, chains, format_flag,
                final_filename, residue_num,
                strand_dist=strand_dist, sheet_dist=sheet_dist, sheet_shift=sheet_shift)
        
    elif classes in [5, 6, 7, 8]:
        packpep_antiparallel(single_pep, length, classes, shift, chains, format_flag,
                             final_filename, residue_num,
                             strand_dist=strand_dist, sheet_dist=sheet_dist, sheet_shift=sheet_shift)


if __name__ == "__main__":
    main()


