#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/BinXBioLab/Data_Science_Capstone_Fall_2022/blob/main/CCI_DataExtraction_25Gb.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# #Cell-Cell Interactions - Notebook 1/2
# ##Notebook purpose: Create input files for CellPhoneDB
# 
# 
# Should version conflicts arise, please see `dataextraction_requirements.txt` for full list of package versions successfully tested.
# 
# 

# ###Setup Notebook for data extraction

# To Do: Click "RESTART RUNTIME" button at end of code output after running pip install cell below!

# In[ ]:


#Need to install scanpy for data import and extraction process, this requires restarting runtime 
get_ipython().system('pip install scanpy==1.9.1')


# In[ ]:


#mount access to google drive to access relevant data for input and output files
from google.colab import drive
drive.mount('/mntDrive')


# TO DO: update directory path below to match your drive and desired working folder.
# 
# Example Syntax:
# 
# 
# 
# ```
# %cd /mntDrive/MyDrive/Capstone/CellphoneDB_Files/22qAGGR01/
# ```
# 
# 

# In[ ]:


#change to relevant directory in mounted drive - you can now assume you are working from this directory
get_ipython().run_line_magic('cd', '/mntDrive/MyDrive/Folder1/Folder2')


# Script tested with the following versions:
# 
# - pandas==1.3.5
# - scanpy==1.9.1

# In[ ]:


#import packages needed for data extraction
import pandas as pd
import scanpy as sc
import os 


# ##Import Data and Transform

# TO DO: fill filepath variable in with your desired input file, .h5ad format expected. 
# 
# Example syntax:
# 
# 
# 
# ```
# adata_filepath = '../pasca.log1p_liger_med_singleR_noglyc_anndata080.h5ad'
# ```
# 
# 

# In[ ]:


adata_filepath = '</folder/folder/data_analysis_output_file.h5ad>' 
adata = sc.read(adata_filepath)
#visual inspection of file format -- may be needed to specify cell type label layer
adata


# Create Counts Data File: 
# 
# columns = cells, rows = genes (ensemble IDs) 

# In[ ]:


#extracts counts layer of anndata file
df_counts = adata.to_df(layer="counts")
#transpose dataframe to fit column & row requirements
df_counts = df_counts.T
# # Set cell ids as column headers
df_counts.columns = adata.obs.index
# # Set rows/genes to Ensemble ID
df_counts.set_index(adata.var.gene_ids, inplace=True) 
# # Visual inspection for correct format
df_counts.head()


# Data cleaning steps to group by patient type: case, control, timepoint, and any other criteria desired.
# 
# TO DO: insert correct obs layers from anndata file that you would like to filter by inside of "< >" 

# In[ ]:


#separate analysis groups, i.e. control vs. case/mutant
#create dataframe with cell, batch (genotype group), and timepoint to split 
df_group = pd.DataFrame(data={'Cell': list(adata.obs['<genotype>'].index), 
                             'Group': list(adata.obs['<genotype>']), 
                             'Timepoint': list(adata.obs['<timepoint>'])})
#visualize dataframe
df_group.head()


# TO DO: Ensure column and filter term match your data & update naming conventions of lists as appropriate. 
# 
# example: `df_group[(df_group['Group'] == 'Control')& (df_group['Timepoint'] == '70d')] ` 
# 
# --> 'Group' & 'Timepoint' are columns and 'Control' & '70d' are the filters in example dataset. 

# In[ ]:


#create lists of cells belonging to each group of interest, in example below we are splitting by control/case and timepoints
control_70d = list(df_group[(df_group['Group'] == 'Control') & (df_group['Timepoint'] == '70d')]['Cell'])
control_150d = list(df_group[(df_group['Group'] == 'Control') & (df_group['Timepoint'] == '150d')]['Cell'])
case_70d = list(df_group[(df_group['Group'] == 'Case') & (df_group['Timepoint'] == '70d')]['Cell'])
case_150d = list(df_group[(df_group['Group'] == 'Case') & (df_group['Timepoint'] == '150d')]['Cell'])


# TO DO: Update naming convention of dataframes as appropriate and update lits within bracket to match naming convention created in previous step!

# In[ ]:


#create table per filter
df_counts_case_70 = df_counts[case_70d]
df_counts_control_70 = df_counts[control_70d]
df_counts_case_150 = df_counts[case_150d]
df_counts_control_150 = df_counts[control_150d]


# Prepare Meta Data File:
# 
# 2 columns: cell & cell type label

# To Do: find correct layer in 'obs' that houses cell to cell type mapping. Leverage adata file organization printed above.

# In[ ]:


adata.obs['<layer that specifies cell type label>']


# In[ ]:


## generate meta file
df_meta = pd.DataFrame(data={'Cell': list(adata.obs['<layer that specifies cell type label>'].index), 
                             'cell_type': list(adata.obs['<layer that specifies cell type label>'])})

df_meta.set_index('Cell',inplace=True)

#check for missing values
print('Missing Values Count:', df_meta.isna().sum().sum())

#visual inspection of table
df_meta.head()


# If Missing Values Count >0, follow next steps, otherwise skip to filtering step "Meta table cells must match..."
# 
# To Do: Ensure shape matches between meta table # cells & counts table # cells 

# In[ ]:


#create list of cells missing celltype label
drop_cells = list(df_meta[df_meta['cell_type'].isna()].index)

#drop cells missing celltype label from meta dataframe
df_meta = df_meta.dropna()
print('Shape of Meta Table:', df_meta.shape)

#update counts table to remove empty cell type
df_counts = df_counts.drop(columns=drop_cells)
print('Shape of Counts Table:', df_counts.shape)


# Meta table cells must match cells in each counts dataset, so we must filter per grouping created above for meta files as well. 
# 
# TO DO: update list inside of `df_meta.index.isin(HERE!!)` to match list naming convention you created in prepping counts data section. 
# 
# 
# 

# In[ ]:


#filter meta table per analysis grouping
df_meta_control_70 = df_meta[df_meta.index.isin(control_70d)]
df_meta_case_70 = df_meta[df_meta.index.isin(case_70d)]
df_meta_control_150 = df_meta[df_meta.index.isin(control_150d)]
df_meta_case_150 = df_meta[df_meta.index.isin(case_150d)]


# ##Save Extracted Data to Google Drive
# 
# To Do:
# Specify output path! These files will be the input for CellPhoneDB. 
# 
# Can take up to 1-1.5 hours to save for large files.

# In[ ]:


output_filepath = '</base_folder/output_folder/>'


# TO DO: If you have created more than 2 groupings, you will need to copy and paste the saving code to also save your other filtered groups as needed. Update the name of the inside the '' within the code to match your filter appropriately but ensure it still contains meta or count for next step usage in cellphonedb!
# 
# The syntax should be as follows:
# - Meta file:
# `<name_of_df_created>.to_csv(os.path.join(output_filepath, '<name_of_new_file_meta_data.txt'), sep='\t') `
# - Counts file:
# <name_of_df_created>.to_csv(os.path.join(output_filepath, '<name_of_new_file_counts_data.h5'), key='counts') 

# In[ ]:


#saves as counts as .h5 file and meta as .txt file, tab delimited - runs within a few minutes
## example of saving case & control for 70d timepoint
df_meta_case_70.to_csv(os.path.join(output_filepath, 'case_70_meta_data.txt'), sep='\t')
df_meta_control_70.to_csv(os.path.join(output_filepath, 'control_70_meta_data.txt'), sep='\t')
df_counts_case_70.to_hdf(os.path.join(output_filepath, 'case_70_counts_data.h5'), key='counts') 
df_counts_control_70.to_hdf(os.path.join(output_filepath, 'control_70_counts_data.h5'), key='counts')

## example of saving case and control for 150d timepoint
df_meta_case_150.to_csv(os.path.join(output_filepath, 'case_150_meta_data.txt'), sep='\t')
df_meta_control_150.to_csv(os.path.join(output_filepath, 'control_150_meta_data.txt'), sep='\t')
df_counts_case_150.to_hdf(os.path.join(output_filepath, 'case_150_counts_data.h5'), key='counts') 
df_counts_control_150.to_hdf(os.path.join(output_filepath, 'control_150_counts_data.h5'), key='counts')

