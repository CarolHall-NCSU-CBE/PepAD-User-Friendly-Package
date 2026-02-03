#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.ticker import MultipleLocator
import argparse


# In[15]:


def rhs(line):
    return line.split("=", 1)[1].strip()

input_signal = True
try:
    sheet_num  = 2
    start = 1
    sheetmove = 0.6
    with open('input.txt', 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.lstrip().startswith("#")]
        
        gnum        = int(rhs(lines[0]))
        pep_num     = int(rhs(lines[1]))
        pep_name    = rhs(lines[2])

        end         = int(rhs(lines[6]))
        ekt_seq     = float(rhs(lines[7]))
        ekt_sheet   = float(rhs(lines[8]))

        lam         = float(rhs(lines[10]))

        hydrophobic = int(rhs(lines[11]))
        polar       = int(rhs(lines[12]))
        charged     = int(rhs(lines[13]))
        other       = int(rhs(lines[14]))
      
        
except Exception as input_err:
    print(f"An error occurred: {input_err}")
    input_signal = False
    


# In[16]:


# assgin headers and read energy profile
headers=["step","Sequence","Score","G_bind","Pagg","rmsd"]
df = pd.read_csv('energyprofile.txt', sep=r'\s+', header=None, names=headers)
df = df.fillna(0)


# In[18]:


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




