#!/usr/bin/env python
# coding: utf-8

# In[16]:


import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.ticker import MultipleLocator
import argparse


# In[17]:


input_signal = True
try:
    with open('input.txt', 'r') as file:
        lines = file.readlines()  # Read all lines into a list

        numbers = lines[1].strip().split()
        line = list(map(float, numbers))
        gnum, pep_num, sheet_num = line

        numbers = lines[3].strip().split()
        pep_name = numbers[0]

        numbers = lines[5].strip().split()
        line = list(map(int, numbers))
        start, end = line

        numbers = lines[7].strip().split()
        line = list(map(float, numbers))
        seed_switch, rand_seed, ekt_seq = line
        
        numbers = lines[9].strip().split()
        line = list(map(float, numbers))
        sheetmove, ekt_sheet, interval = line # sheet move on > 0, off = 0
        
        numbers = lines[11].strip().split()
        line = list(map(float, numbers))
        rmsdx, rmsdy, rmsdz, dx, dy, dz = line # sheet move conditions

        numbers = lines[13].strip().split()
        line = list(map(float, numbers))
        lam = line[0]

        numbers = lines[15].strip().split()
        line = list(map(int, numbers))
        hydrophobic, polar, charged, other = line
        
except Exception as input_err:
    print(f"An error occurred: {input_err}")
    input_signal = False


# In[18]:


# assgin headers and read energy profile
headers=["step","Sequence","Score","E_bind","S_bind","I_hydrophobic","I_propensity"]
df = pd.read_csv('energyprofile.txt', sep=r'\s+', header=None, names=headers)
df_rmsd = pd.read_csv('rmsd.txt', sep=r'\s+')

# calculate free energy and Pagg
df["G_bind"]=df["E_bind"]-df["S_bind"]
df["Pagg"]=df["I_hydrophobic"]+df["I_propensity"]

# Merge the two dataframes on the 'step' column
df = pd.merge(df, df_rmsd, on='step', how='outer')

# Fill NaN values with zero
df = df.fillna(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PepAD peptide analysis")
    parser.add_argument("--top", type=int, default=10, help="Number of top peptides to save")
    parser.add_argument("--plot", choices=["score", "rmsd", "both", "none"], default="both", help="select: score, rmsd, both, none")
    parser.add_argument("--score_rolling",type=int, default=0, help="score rolling average range, need int > 0")
    parser.add_argument("--rmsd_rolling", type=int, default=0, help="rmsd rolling average range, need int > 0")
        
    args = parser.parse_args()

    best_unique_pep = args.top
    plot_type       = args.plot
    score_r         = args.score_rolling
    rmsd_r          = args.rmsd_rolling

    colors = ["#DB594B", "#4B89DB", "#DB834B", "#4BD5DB", "#866C5A", "#614841", "#854BDB"]
    if (plot_type=="score"):
        fig, axs = plt.subplots(1, figsize=(10,6), constrained_layout=True)
        
        
        if (score_r > 1):
#             df['Rolling'] = df['Score'].rolling(score_r).mean()
            axs.plot(df['step'], df['Score'].rolling(score_r).mean(), label = "Score (kcal/mol)", color = colors[0], linewidth=1 )
        else:
            axs.plot(df['step'], df["Score"], label = "Score (kcal/mol)", color = colors[0], linewidth=1 )
            
        axs.set_xlabel('Step', fontsize=18)
        axs.set_ylabel('{} (kcal/mol)'.format("Score"), fontsize=20)
        axs.legend(loc="lower center",ncol=2,
                   bbox_to_anchor=(0.5, 1.01),fontsize=20,frameon=False )
        axs.tick_params(axis='both', labelsize=18)
        axs.set_xlim(0,end)
        
    elif (plot_type=="rmsd"):
        fig, axs = plt.subplots(1, figsize=(10,6), constrained_layout=True)
        
        if (rmsd_r > 1):
#             df['Rolling'] = df['rmsd'].rolling(rmsd_r).mean()
            axs.plot(df['step'], df['rmsd'].rolling(rmsd_r).mean(), label = "RMSD (\u212B)", color = colors[1], linewidth=1 )
        else:
            axs.plot(df['step'], df["rmsd"], label = "RMSD (\u212B)", color = colors[1], linewidth=1 )
            
        axs.set_xlabel('Step', fontsize=18)
        axs.set_ylabel('RMSD (\u212B)', fontsize=20)
        axs.legend(loc="lower center",ncol=2,
                   bbox_to_anchor=(0.5, 1.01),fontsize=20,frameon=False )
        axs.tick_params(axis='both', labelsize=18)
        axs.set_xlim(0,end)
    
    elif (plot_type=="both"):
        fig, axs = plt.subplots(1, figsize=(10,6), constrained_layout=True)
        
        if (score_r > 1):
#             df['Rolling_score'] = df['Score'].rolling(score_r).mean()
            axs.plot(df['step'], df['Score'].rolling(score_r).mean(), label = "{} (kcal/mol)".format("Score"), color = colors[0], linewidth=1 )        
        else:
            axs.plot(df['step'], df["Score"], label = "Score (kJ/mol)", color = colors[0], linewidth=1 )
        
        axs.set_xlabel('Step', fontsize=18)
        axs.set_ylabel('Score (kcal/mol)', fontsize=20)
        
        ax2 = axs.twinx()
        if (score_r > 1):
#             df['Rolling_rmsd'] = df['rmsd'].rolling(rmsd_r).mean()
            ax2.plot(df['step'], df['rmsd'].rolling(rmsd_r).mean(), label = "RMSD (\u212B)", color = colors[1], linewidth=1 )
        else:
            ax2.plot(df['step'], df["rmsd"], label = "RMSD (\u212B)", color = colors[1], linewidth=1 )
        
        ax2.set_ylabel('RMSD (\u212B)', fontsize=20)
        
        lines_1, labels_1 = axs.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        axs.legend(lines_1 + lines_2, labels_1 + labels_2,loc="lower center",ncol=2,
                   bbox_to_anchor=(0.5, 1.01),fontsize=20,frameon=False )

        ax2.set_ylabel('RMSD (\u212B)', fontsize=20)
        axs.set_xlim(0,end)

        axs.tick_params(axis='both', labelsize=18)
        ax2.tick_params(axis='both', labelsize=18)
    
    plt.savefig('step_evolution', dpi=600)


# In[19]:


# Find the index of the row with the minimum value in the column
score_min_index = df['Score'].idxmin()

df_sorted = df.sort_values(by='Score', ascending=True) # ascending=True from small to large, =false from large to small
pep_counts = df_sorted['Sequence'].value_counts() # counting duplicates
pep_counts = pd.DataFrame(pep_counts)

# Score of each unique sequence is averaged
df_unique_ave = df.groupby("Sequence").mean() #根据sequence，每个相同sequence求平均
df_unique_ave = df_unique_ave.sort_values(by='Score', ascending=True)
df_unique_ave = df_unique_ave.merge(pep_counts, how='left', left_on='Sequence', right_index=True)
df_unique_ave = df_unique_ave.rename(columns={'count': 'Counts'})
df_unique_ave = df_unique_ave.drop('step', axis=1) # axis = 0 is row, = 1 is column
df_unique_ave = df_unique_ave.reset_index()
df_ave_unique = df_unique_ave.head(best_unique_pep)

# Drop duplicate sequences, keeping the row with the lowest score for each sequence
df_unique_best = df_sorted.drop_duplicates(subset='Sequence', keep='first')

# powerful merge, how="left" is left join method, all df_unique is included, the pep_counts is matched,
# left_on='Sequence' is the left matching column, right_index=True is the index of pep_counts used for matching
df_unique_best = df_unique_best.merge(pep_counts, how='left', left_on='Sequence', right_index=True)
df_unique_best = df_unique_best.rename(columns={'count': 'Counts'})

# Select the top 5 rows with the lowest scores among the unique sequences
df_min_unique = df_unique_best.head(best_unique_pep)


# In[20]:


# write out the result
filename="PepAD report.txt"
mean_values = df.mean(numeric_only=True).to_frame().T # transform to data frame, and transpose 
std_values = df.std(numeric_only=True).to_frame().T
mean_values.index = ['Ave'] # rename the index
std_values.index = ['Std']
mean_std_df = pd.concat([mean_values, std_values]) # combine Ave and Std column
mean_std_df = mean_std_df.drop(columns=['step']) # remove "Step" column

with open ( filename,'w') as f:
############### Properties #####################
    f.write("---{} unique peptides with best score (energy profile)---\n".format(best_unique_pep))
    f.write(df_min_unique.to_string(index=False, float_format='%.2f'))
    f.write("\n\n")
    





# In[ ]:





# In[ ]:




