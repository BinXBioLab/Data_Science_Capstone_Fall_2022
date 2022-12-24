#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/BinXBioLab/Data_Science_Capstone_Fall_2022/blob/main/CCI_Run_CellPhoneDB_25Gb.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# #Cell-Cell Interactions - Notebook 2/2
# ##Notebook purpose: Run CellPhoneDB Statistical Analysis Method
# 
# Should version conflicts arise, please see `run_cpdb_requirements.txt` for full list of successfully tested package versions. 
# 
# 

# ###Setup Notebook for CellphoneDB

# pip install cellphonedb will ensure all dependency packages match requirements set by CellPhoneDB package. 
# 
# Successfully run on version 3.1.0 - suggested to use this version.
# 
# **TO DO:**
# 
# Restart runtime required after pip install of cellphonedb!

# In[ ]:


get_ipython().system('pip install cellphonedb==3.1.0')


# In[ ]:


import cellphonedb
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
#helpful to navigate to folder where input files are located
get_ipython().run_line_magic('cd', '/mntDrive/MyDrive/Folder1/Folder2/Samples')


# ##Run Statistical Analysis

# TO DO: Create output folders for each analysis you are running
# 
# Example below for analysis run across control and case at two time points: 
# 
# 
# 
# ```
# !mkdir CONTROL70_out
# !mkdir CONTROL150_out
# !mkdir CASE70_out
# !mkdir CASE150_out
# ```
# 
# 

# In[ ]:


get_ipython().system('mkdir Output_folder1')
get_ipython().system('mkdir Output_folder2')


# TO DO: locate input data file paths with respect to current directory (meta_data.txt & counts_data.txt)
# 
# Expected Output: Verbose outputs will print periodically. Statistical Analysis output files will be saved in folder you specify in the `--output-path=` parameter.
# 
# Note: Can take up to 3.5 hours to run. If verbose output ends with "^C" memory limit has been reached and larger RAM space is required for input data files.
# 
# Syntax Example on previous run:
# 
# 
# ```
# !cellphonedb method statistical_analysis control_70_meta_data.txt control_70_counts_data.h5 --output-path=/mntDrive/MyDrive/Capstone/CellphoneDB_Files/22qAGGR01/CONTROL70_out
# ```
# 
# 

# In[ ]:


get_ipython().system('cellphonedb method statistical_analysis group1_meta_data.txt group1_counts_data.h5 --output-path=/TestCode')


# In[ ]:


get_ipython().system('cellphonedb method statistical_analysis group2_meta_data.txt group2_counts_data.h5 --output-path=/Output_folder2')


# ##Generate Visualizations - Dot Plot & Heatmap Plot

# To Do: Ensure proper path is shown below to output files
# 
# Note:
# To limit size of dotplot visualization include files with the following arguments: 
# 
# --rows: File with a list of rows to plot, one per line (defualt will use all available)
# 
# --columns: File with a list of columns to plot, one per line (default will use all available)
# 
# ----------------------
# 
# --output-name is used to save file as jpg (default is pdf, remove argument if you wish to save as pdf)
# 
# 
# Example of previously run command:
# 
# 
# ```
# ##plot case 
# !cellphonedb plot dot_plot --means-path CASE_out/means.txt --pvalues-path CASE_out/pvalues.txt --output-path=/mntDrive/MyDrive/Capstone/CellphoneDB_Files/22qAGGR01/CASE_out
# ```
# 
# 
# 
# 

# In[ ]:


#plot dotplot group 1
get_ipython().system('cellphonedb plot dot_plot --output-name dotplot.jpg --means-path Output_folder1/means.txt --pvalues-path Output_folder1/pvalues.txt')


# In[ ]:


#plot dotplot group 2
get_ipython().system('cellphonedb plot dot_plot --output-name dotplot.jpg --means-path Output_folder2/means.txt --pvalues-path Output_folder2/pvalues.txt')


# Heatmap plot requires R package installation, only run if heatmap plotting is desired.

# In[ ]:


#installing "pheatmap" R package requires the following version changes to be made to python packages
##upgrade to pandas==1.3.5
##downgrade to rpy2==3.5.1
#Restart Runtime required after version changes are made!

get_ipython().system('pip install pandas==1.3.5')
get_ipython().system('pip install rpy2==3.5.1')

import rpy2


# In[ ]:


#reset current working directory after runtime is restarted -- check to ensure you are pointing to folder where input files and "out/" folder are located
get_ipython().run_line_magic('cd', '/mntDrive/MyDrive/Capstone/CellphoneDB_Files/')


# In[ ]:


#import r2py for notebook
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[ ]:


#import R package for heatmap
get_ipython().run_line_magic('R', "install.packages('pheatmap')")


# To Do: Ensure proper path for meta_data (generated in previous script) and pvalues.txt file (generated by cellphonedb)
# 
# Note: --count-name ensures plot is saved as jpg (defualt is pdf, remove argument if you wish to save as pdf)

# In[ ]:


#heatmap for group 1
get_ipython().system('cellphonedb plot heatmap_plot meta_data.txt --pvalues-path Output_folder1/pvalues.txt')


# In[ ]:


#heatmap for group 2
get_ipython().system('cellphonedb plot heatmap_plot meta_data.txt --pvalues-path Output_folder2/pvalues.txt')

