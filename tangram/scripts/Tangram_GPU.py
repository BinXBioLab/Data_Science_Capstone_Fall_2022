#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from google.colab import drive
drive.mount('/content/drive')
get_ipython().system('mkdir /content/MyDrive')
get_ipython().system('mount --bind /content/drive/My\\ Drive /content/MyDrive')
# the magic % will help with the cd command :)
get_ipython().run_line_magic('cd', '/content')


# In[ ]:


get_ipython().system('pip install scanpy')
get_ipython().system('pip install squidpy')
get_ipython().system('pip install scikit-image')
get_ipython().system('pip install tangram-sc')
get_ipython().system('pip install gdown')
# to upgrade
get_ipython().system('pip install --upgrade gdown')


# In[12]:


import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg


# In[2]:


get_ipython().run_line_magic('cd', "'/content/MyDrive/spatial_scRNAseq_test'")


# In[3]:


#!gdown --fuzzy 'https://drive.google.com/file/d/19O5HyRPUcaeYZut7m_wC7ND6Tb_vNVle/view?usp=share_link'


# # 1. Load in Files

# In[13]:


adata_sc = sc.read_h5ad("pasca.log1p_liger_med_singleR_noglyc.h5ad")


# In[57]:


adata_st = sq.datasets.visium('V1_Human_Brain_Section_1')   
adata_st = sq.datasets.sc_mouse_cortex()
#adata_st = sq.datasets.visium('V1_Human_Brain_Section_2') 


# In[15]:


img = sq.datasets.visium_fluo_image_crop()
#img = sq.datasets.visium_hne_image_crop()


# # 2. Preprocess spatial data

# In[58]:


sc.pp.calculate_qc_metrics(adata_st, inplace=True)


# In[59]:


adata_st


# In[60]:


sc.pp.normalize_total(adata_st, inplace=True)
sc.pp.log1p(adata_st)
sc.pp.highly_variable_genes(adata_st, flavor="seurat", n_top_genes=2000, inplace=True)


# In[61]:


sc.pp.neighbors(adata_st)
sc.tl.umap(adata_st)
sc.tl.leiden(adata_st, key_added="clusters")


# In[47]:


# Pre processing of adata_st ends


# In[48]:


#adata_st.write_h5ad("adata_st_postpre.h5ad")


# In[49]:


#adata_st = sc.read_h5ad("adata_st_postpre.h5ad")


# # 3. Generate Map

# ## 3.1 Select Genes

# In[62]:


adata_sc.uns['log1p']["base"] = None


# ### 3.1.1 Rank Gene Groups (Recommended)

# In[63]:


sc.tl.rank_genes_groups(adata_sc, groupby="nowakowski_med", use_raw=False)
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
len(markers)


# ### 3.1.2 Highly variable Genes (Alternative)

# In[25]:


'''
sc.pp.highly_variable_genes(adata_sc)
highvar_df = pd.DataFrame(adata_sc.var["highly_variable"])
highvar_df[highvar_df["highly_variable"] == True]
highvar_markers = list(highvar_df[highvar_df["highly_variable"] == True].index.values)
print(highvar_markers[:20])
'''


# ## 3.2 Single Anndata

# ### 3.2.1 Single Preprocess AnnData

# In[9]:


tg.pp_adatas(adata_sc, adata_st, genes=markers)


# ### 3.2.2 Single Generate Map

# This is the most important process

# In[11]:


import torch
torch.cuda.empty_cache()
ad_map = tg.map_cells_to_space(adata_sc, adata_st,
    mode="cells",
#     mode="clusters",
#     cluster_label='cell_subclass',  # .obs field w cell types
    density_prior='rna_count_based',
    num_epochs=250,
#     device="cuda:0",
    device='cpu',
)


# In[23]:


#adata_st
#adata_st.write_h5ad("adata_st_postmap.h5ad")


# In[ ]:


#adata_st = sc.read_h5ad("adata_st_postmap.h5ad")
#ad_map = sc.read_h5ad("ad_map.h5ad")


# ### 3.2.3 Save Map

# In[26]:


adata_st.write_h5ad("adata_st.h5ad")
ad_map.write_h5ad("ad_map.h5ad")


# ## 3.3 Double Anndata

# ###  3.3.1 Separate to two single cell data

# In[64]:


adata_sc.obs['sampleID']


# In[65]:


filter_val = ['GSM4306931_Control_1', 'GSM4306932_Control_2']
adata_sc1 = adata_sc[adata_sc.obs['sampleID'].isin(filter_val),:]
adata_sc2 = adata_sc[~adata_sc.obs['sampleID'].isin(filter_val),:]


# ### 3.3.2 Double Preprocess Data

# In[66]:


tg.pp_adatas(adata_sc1, adata_st, genes=markers)
tg.pp_adatas(adata_sc2, adata_st, genes=markers)


# In[56]:


adata_sc2


# ### 3.3.3 Double Generate Map

# In[68]:


epochs = 500

import torch
torch.cuda.empty_cache()
ad_map1 = tg.map_cells_to_space(adata_sc1, adata_st,
    mode="cells",
#     mode="clusters",
#     cluster_label='cell_subclass',  # .obs field w cell types
    density_prior='rna_count_based',
    num_epochs= epochs,
     device="cuda:0",
#    device='cpu',
)
torch.cuda.empty_cache()
ad_map2 = tg.map_cells_to_space(adata_sc2, adata_st,
    mode="cells",
#     mode="clusters",
#     cluster_label='cell_subclass',  # .obs field w cell types
    density_prior='rna_count_based',
    num_epochs= epochs,
     device="cuda:0",
#    device='cpu',
)


# ### 3.3.4 Double Save Map

# In[14]:


adata_st


# In[15]:


adata_st.write_h5ad("adata_st.h5ad")
ad_map1.write_h5ad("ad_map1.h5ad")
ad_map2.write_h5ad("ad_map2.h5ad")


# In[ ]:




