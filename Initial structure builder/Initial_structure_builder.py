#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
from PeptideBuilder import Geometry
import PeptideBuilder
import pandas as pd
import scipy.linalg
import matplotlib.pyplot as plt
import argparse
import Bio.PDB
import hydride
import subprocess


# In[8]:


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

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)
    
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
    
    
    df = df[(df['atom_name'] != "HXT")]                      # remove H
    df = df[(df['atom_name'] != "H2") | (df['resid'] != 1)]  # remove H2
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

        HA= df[(df['atom_name'] == 'HA') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
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

        HA= df[(df['atom_name'] == 'HA') & (df['resid'] == last_resid)].iloc[0][['x', 'y', 'z']].to_numpy()
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
    
    v1 = N - C; V2 = H - N
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



# In[9]:


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
    bond_length = length of N–H bond
    bond_angle_deg = CA–N–H angle (in degrees)
    dihedral_deg = CB–CA–N–H dihedral angle (in degrees)
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


# In[10]:


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


# In[33]:


def packpep(single_peptite, chain_length, classes, shifts, chains_per_sheet, format_flag, o_filename, 
            core_residues="e", strand_dist = 4.8, sheet_dist = 11.5, sheet_shift = 0):
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

 


# In[35]:


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
        print("Cannot determine Class. Stop.")
        exit()

    ###########################################
    ########## create single pep ##############
    ###########################################
    temp_sequence = sequence
    geo = Geometry.geometry(temp_sequence[0])      
    structure = PeptideBuilder.initialize_res(geo)   # create the first residue
    
    for aa in temp_sequence[1:]:                     # create the rest residues
        geo = Geometry.geometry(aa)
        geo.phi = angles[0]
        geo.psi_im1 = angles[1]
        PeptideBuilder.add_residue(structure, geo)
    # add terminal oxygen (OXT) to the final glycine
    PeptideBuilder.add_terminal_OXT(structure)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save("temp_sequence.pdb")

    ###########################################
    ############# add hydrogen ################
    ###########################################
    input_file = "temp_sequence.pdb"  
    output_file = "output_structure.pdb" 
    
    # Run the hydride (Add hydrogen)
    try:
        subprocess.run(
            ["hydride", "-i", input_file, "-o", output_file],
            check=True  # Ensures it raises an error if the command fails
        )
        print(f"Hydrogens added successfully. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running hydride: {e}")

    ###########################################
    ############# Modify termini ##############
    ###########################################
    if cap_flag == 0:      # Add hydrogen at N-terminus and remove H at C-terminus
        add_NH_remove_OH("output_structure.pdb", "output_structure2.pdb")
        length=len(sequence)
        
    elif cap_flag in [1, 2]:    # Add Caps to peptide
        add_caps("output_structure.pdb","output_structure2.pdb", cap_flag)
        length=len(sequence) + 2
    else:
        print("Cannot determine caps or not. Stop")
        exit()
        
    ###########################################
    ############# PEP ALIGNMENT ###############
    ###########################################
    single_pep = alignment(length,1,"output_structure2.pdb", "output_aligned.pdb", order = 1)
    
    ###########################################
    ############## Pack sheets ################
    ###########################################
    if classes in [1, 2, 3, 4]: #parallel sheets
        packpep(single_pep, length, classes, shift, chains, format_flag,
                final_filename, residue_num, strand_dist, sheet_dist, sheet_shift)
        print(f"Packed {sequence} into {final_filename}")
    elif classes in [5, 6, 7, 8]:
        packpep_antiparallel(single_pep, length, classes, shift, chains, format_flag,
                final_filename, residue_num, strand_dist, sheet_dist, sheet_shift)
        print(f"Packed {sequence} into {final_filename}")


if __name__ == "__main__":
    main()


# In[ ]:


def main2(sequence, pep_name, residue_num):
#     parser = argparse.ArgumentParser(description="Peptide builder")

#     parser.add_argument('-seq', type=str, required=True, help="Peptide sequence (e.g., GNNQQNY)")
#     parser.add_argument('-c', type=int, required=True, choices=range(1, 9), help="Class (1-8 only)")
#     parser.add_argument('-sh', type=float, default=0.0, help="Shift value: -1, 0, or 1")
#     parser.add_argument('-n', type=int, default=8, help="Number of chains (default: 8)")
#     parser.add_argument('-p', type=int, default=0, help="terminal patch uncap = 0 (default), capped=1")
#     parser.add_argument('-r', type=str, default="e", help="residue number to be packed in the sheets even = e default, odd = o")
#     args = parser.parse_args()

#     print("Sequence:", args.seq)
#     print("Class:", args.c)
#     print("Shift:", args.sh)
#     print("Number of chains:", args.n)
#     print("Cap flag:", args.p)
   
    # sequence = "NNNGTLNIAITVNN"
    classes  = 1
    shift    = -1
    chains   = 12
    cap_flag = 2
    # residue_num = "e"
    format_flag = 1
    output_name = pep_name + ".pdb"
    
    
    
    
    if classes in [1, 2, 3, 4]: #parallel sheets
        angles = [-119, 113]
    elif classes in [5, 6, 7, 8]:
        angles = [-139, 135]
    else:
        print("Cannot determine Class. Stop.")
        exit()

    ###########################################
    ########## create single pep ##############
    ###########################################
    temp_sequence = sequence
    geo = Geometry.geometry(temp_sequence[0])      
    structure = PeptideBuilder.initialize_res(geo)   # create the first residue
    
    for aa in temp_sequence[1:]:                     # create the rest residues
        geo = Geometry.geometry(aa)
        geo.phi = angles[0]
        geo.psi_im1 = angles[1]
        PeptideBuilder.add_residue(structure, geo)
    # add terminal oxygen (OXT) to the final glycine
    PeptideBuilder.add_terminal_OXT(structure)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save("temp_sequence.pdb")

    ###########################################
    ############# add hydrogen ################
    ###########################################
    input_file = "temp_sequence.pdb"  
    output_file = "output_structure.pdb" 
    
    # Run the hydride (Add hydrogen)
    try:
        subprocess.run(
            ["hydride", "-i", input_file, "-o", output_file],
            check=True  # Ensures it raises an error if the command fails
        )
        print(f"Hydrogens added successfully. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running hydride: {e}")

    ###########################################
    ############# Modify termini ##############
    ###########################################
    if cap_flag == 0:      # Add hydrogen at N-terminus and remove H at C-terminus
        add_NH_remove_OH("output_structure.pdb", "output_structure2.pdb")
        length=len(sequence)
        
    elif cap_flag in [1, 2]:    # Add Caps to peptide
        add_caps("output_structure.pdb","output_structure2.pdb", cap_flag)
        length=len(sequence) + 2
    else:
        print("Cannot determine caps or not. Stop")
        exit()
        
    ###########################################
    ############# PEP ALIGNMENT ###############
    ###########################################
    single_pep = alignment(length,1,"output_structure2.pdb", "output_aligned.pdb", order = 1)
    
    ###########################################
    ############## Pack sheets ################
    ###########################################
    if classes in [1, 2, 3, 4]: #parallel sheets
        packpep(single_pep, length, classes, shift, chains, format_flag,
                output_name, residue_num)
        
    elif classes in [5, 6, 7, 8]:
        packpep_antiparallel(single_pep, length, classes, shift, chains, format_flag,
                             output_name, residue_num)






