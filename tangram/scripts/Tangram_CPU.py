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


# In[1]:


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


# # 1. Load in Files

# There are 3 files you have to load in to this notebook.
# 
# 1. `adata_sc`: The single cell data file
# 
# 2. `adata_st`: The preprocessed spatial data file. It is a result of the previous notebook `Tangram_GPU.ipynb`. 
# 
# 3. `ad_map`: The spatial mapping file. It is a result of the previous notebook `Tangram_GPU.ipynb`.

# In[2]:


adata_sc = sc.read_h5ad("pasca.log1p_liger_med_singleR_noglyc.h5ad")
adata_st = sc.read_h5ad("adata_st_postproj.h5ad")
ad_map = sc.read_h5ad("ad_map.h5ad")
#ad_ge = sc.read_h5ad("ad_ge.h5ad")


# # 2. Analyzing data

# ## 2.1 Visualizing initial data 

# In[3]:


img = sq.datasets.visium_fluo_image_crop()


# In[4]:


sc.pl.umap(
    adata_st, color=["clusters"], palette=sc.pl.palettes.default_20
)


# In[5]:


clusters_colors = dict(
    zip([str(i) for i in range(18)], adata_st.uns["clusters_colors"])
)


# In[6]:


plt.rcParams["figure.figsize"] = (8, 8)
ad = adata_st.copy()
sc.pl.spatial(
  ad,
  img_key="hires",
  color=["clusters","n_genes_by_counts"],
  size=1.5,
  palette=[
    v
    for k, v in clusters_colors.items()
    if k in ad.obs.clusters.unique().tolist()
    ],
  legend_loc='right margin',
  show=False,
    )

plt.tight_layout()


# In[7]:


fig, axs = plt.subplots(1, 2, figsize=(20, 5))
sc.pl.spatial(
    adata_st, color="clusters", alpha=0.7, frameon=False, show=False, ax=axs[0]
)
sc.pl.umap(
    adata_sc, color="nowakowski_med", size=10, frameon=False, show=False, ax=axs[1]
)
plt.tight_layout()


# ## 2.2 Pre-process AnnDatas

# This preprocessing is required for `adata_sc` to be in the correct format for further experiments.

# In[8]:


adata_sc.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata_sc, groupby="nowakowski_med", use_raw=False)
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
len(markers)
tg.pp_adatas(adata_sc, adata_st, genes=markers)


# In[9]:


# alternative 2: saving the experiment result from the previous gpu notebook can be used but pp_adatas don't take much time.
# It is not worth generating another intermediate file.
#adata_sc = sc.read_h5ad("adata_sc.h5ad")


# ## 2.3 Visualize cell type maps

# In[10]:


tg.project_cell_annotations(ad_map, adata_st, annotation="nowakowski_med")
annotation_list = list(pd.unique(adata_sc.obs['nowakowski_med']))
tg.plot_cell_annotation_sc(adata_st, annotation_list)


# In[11]:


#separation of two samples visualize.


# In[12]:


tg.plot_training_scores(ad_map, bins=20, alpha=.5)


# ## 2.4 Generate new spatial data via aligned single cells

# In[13]:


# This part to go to CPU.
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
#ad_ge = sc.read_h5ad("ad_ge.h5ad")


# In[14]:


import gc
gc.collect()


# In[15]:


df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_st, adata_sc)
df_all_genes


# In[16]:


tg.plot_auc(df_all_genes)


# ### Inspect predictions

# In[17]:


# we have to use low letter
GOI=['pax6','sox2','foxg1','neurod2','prodh','ranbp1','crkl','hira']


# In[ ]:


tg.plot_genes_sc(GOI, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)


# In[49]:


adata_st


# two kinds of data, one is the patient the other is control and if we can predicted(control) predicted(patient).

# In[18]:


plt.rcParams["figure.figsize"] = (20, 20)
ad = ad_ge.copy()
sc.pl.spatial(
  ad,
  img_key="hires",
  color=GOI,
  size=1.5,
  palette=[
    v
    for k, v in clusters_colors.items()
    if k in ad.obs.clusters.unique().tolist()
    ],
  legend_loc='right margin',
  show=False,
    )

plt.tight_layout()


# ## 2.? Try to compare two groups

# In[19]:


adata_st


# In[20]:


adata_sc


# In[21]:


ad_map


# In[22]:


GOI=['pax6','sox2','foxg1','neurod2','prodh','ranbp1','crkl','hira']
sc.pl.spatial(
            ad_ge,
            color=GOI
        )


# ## 2.5 Deconvolution via alignment

# In[23]:


sf = adata_st.uns['spatial']['V1_Human_Brain_Section_1']["scalefactors"]["tissue_hires_scalef"]
img = sq.im.ImageContainer(adata_st.uns['spatial']['V1_Human_Brain_Section_1']['images']['hires'],
                           scale = sf, library_id = 'V1_Human_Brain_Section_1')


# In[24]:


sq.im.process(img=img, layer="image", method="smooth")
sq.im.segment(                                                                                
    img=img,
    layer="image_smooth",
    method="watershed",
    channel=0,
)


# In[25]:


inset_y = 1500
inset_x = 1700
inset_sy = 400
inset_sx = 500

fig, axs = plt.subplots(1, 3, figsize=(30, 10))
# spatial how?
sc.pl.spatial(
    adata_st, 
    color="rna_count_based_density",   
    alpha=0.7, frameon=False, show=False, ax=axs[0], title=""
)
axs[0].set_title("Clusters", fontdict={"fontsize": 20})
rect = mpl.patches.Rectangle(
    (inset_y * sf, inset_x * sf),
    width=inset_sx * sf,
    height=inset_sy * sf,
    ec="yellow",
    lw=4,
    fill=False,
)
axs[0].add_patch(rect)

axs[0].axes.xaxis.label.set_visible(False)
axs[0].axes.yaxis.label.set_visible(False)

axs[1].imshow(
    img["image"][inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx, 0, 0]
    / 65536,
    interpolation="none",
)
axs[1].grid(False)
axs[1].set_xticks([])
axs[1].set_yticks([])
axs[1].set_title("DAPI", fontdict={"fontsize": 20})

crop = img["segmented_watershed"][
    inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx
].values.squeeze(-1)
crop = skimage.segmentation.relabel_sequential(crop)[0]
cmap = plt.cm.plasma
cmap.set_under(color="black")
axs[2].imshow(crop, interpolation="none", cmap=cmap, vmin=0.001)
axs[2].grid(False)
axs[2].set_xticks([])
axs[2].set_yticks([])
axs[2].set_title("Nucleous segmentation", fontdict={"fontsize": 20});


# In[26]:


features_kwargs = {
    "segmentation": {
        "label_layer": "segmented_watershed",
        "props": ["label", "centroid"],
        "channels": [1, 2],
    }
}
# calculate segmentation features
sq.im.calculate_image_features(
    adata_st,
    img,
    layer="image",
    key_added="image_features",
    n_jobs=4,
    features_kwargs=features_kwargs,
    features="segmentation",
    mask_circle=True,
)


# In[27]:


adata_st.obs["cell_count"] = adata_st.obsm["image_features"]["segmentation_label"]
sc.pl.spatial(adata_st, color=["cell_count"], frameon=False)


# In[28]:


adata_st.obsm['image_features']


# In[34]:


ad = adata_st.copy()
print(ad.var_names)
ad.obs.columns
#sq.pl.spatial_scatter(ad, color=["nowakowski_med"])


# In[35]:


ad


# In[37]:


tg.create_segment_cell_df(ad)
ad.uns["tangram_cell_segmentation"].head()


# In[38]:


ad.obs = pd.concat([adata_st.obs, adata_st.obsm["tangram_ct_pred"]], axis=1)
annotation_list = ad_map.obs['nowakowski_med']
annotation_list


# In[44]:


ad


# In[45]:


ad_map


# In[46]:


adata_sc


# In[75]:


xs = ad.obsm["spatial"][:, 1]
ys = ad.obsm["spatial"][:, 0]
cell_count = ad.obsm["image_features"]["segmentation_label"]

df_segmentation = ad.uns["tangram_cell_segmentation"]
centroids = ad.obsm["tangram_spot_centroids"]

    # create a dataframe
df_vox_cells = df_vox_cells = pd.DataFrame(
    data={"x": xs, "y": ys, "cell_n": cell_count, "centroids": centroids},
    index=list(ad.obs.index),
)
df_vox_cells


# In[60]:


def one_hot_encoding(l, keep_aggregate=False):
    df_enriched = pd.DataFrame({"cl": l})
    for i in l.unique():
        df_enriched[i] = list(map(int, df_enriched["cl"] == i))
    if not keep_aggregate:
        del df_enriched["cl"]
    return df_enriched


# In[95]:


resulting_voxels = np.argmax(ad_map.X, axis=1)
if "F_out" in ad_map.obs.keys():
    filtered_voxels_to_types = [
        (j, adata_sc.obs['nowakowski_med'][k])
        for i, j, k in zip(
            ad_map.obs["F_out"], resulting_voxels, range(len(adata_sc))
        )
        if i > threshold
    ]

    vox_ct = filtered_voxels_to_types

else:
    vox_ct = [(resulting_voxels, adata_sc.obs['nowakowski_med'])]
print(len(resulting_voxels))
print(len(adata_sc.obs['nowakowski_med']))
df_classes = one_hot_encoding(adata_sc.obs['nowakowski_med'])
for index, i in enumerate(df_classes.columns):
    df_vox_cells[i] = 0
colnames = list(df_vox_cells.columns)
print(colnames)
nan_ind = colnames.index(np.nan)
print(ind)
for i in range(0, len(resulting_voxels)):
    k = resulting_voxels[i]
    v = adata_sc.obs["nowakowski_med"][i]
    if(v is not np.nan):
        ind = df_vox_cells.columns.get_loc(v)
    else:
        ind = nan_ind
    df_vox_cells.iloc[k, ind] += 1

ad.obsm["tangram_ct_count"] = df_vox_cells


# In[ ]:


tg.count_cell_annotations(
    adata_map = ad_map,
    adata_sc = adata_sc,
    adata_sp = ad,
    annotation = 'nowakowski_med'
)
# nowakowski_med
adata_st.obsm["tangram_ct_count"].head()


# In[96]:


adata_segment = tg.deconvolve_cell_annotations(ad)
adata_segment.obs.head()


# 
