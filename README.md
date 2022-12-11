# Data Science Capstone Fall 2022
## Authors: Ju, Katharina, Darvesh

# Setting Up
Most files that will be needed are contained within the repository. The only files that are not contained within the repository are those files in `raw_data/Pasca_scRNAseq/inter/`. They are not included because they are too large to commit to GitHub via normal means and we lack adequate storage space and bandwidth on git lfs (Git Large-file Storage). However, these files are also generated via running the pipeline so they are not strictly necessary for running the code from scratch. A `git clone` should be sufficient to get started. If you have access to the Google Drive folder with the data that belongs in `raw_data/Pasca_scRNAseq/inter/`, go ahead and copy it there. The files in there will automatically be ignored via the `.gitignore` file so any changes to them will not be tracked. This may change in the future but it is unlikely due to the extremely large nature of the files (>16 Gb).

## Python
This repo uses python 3.10.4 for all the python code (scripts and notebooks). We recommend creating a conda environment as follows:

```
conda env create -f environment.yml
```

This will install all the necessary python dependencies on your machine. To access this environment, type the following (`scrna` is the default name of the conda environment):

```
conda activate scrna
```

## R
We use version 4.2.1 for R. The R packages are installed in the scripts so there's no need to worry about installing the R packages. This may change later. Please download the appropriate version of R [here](https://archive.linux.duke.edu/cran/) based on your operating system. 

## Docker
Coming soon...



# SpatialTranscriptomics

## General structure

The code is separated into two pipelines: **Tangram_GPU** and **Tangram_CPU**. **Tangram_GPU** will be generating the map that can be used in **Tangram_CPU**. **Tangram_CPU** will be running all the experiments. The only requirement is to choose the same marker genes. 

## Tangram_GPU
### Input
Tangram_GPU only needs one input file: 
- adata_sc : pasca.log1p_liger_med_singleR_noglyc.h5ad.

### Process
In Tangram_GPU, the user has to select which 10xGenomics dataset they will utilize. Choices of visium dataset can be checked in [squidpy](https://squidpy.readthedocs.io/en/stable/api/squidpy.datasets.visium.html#squidpy.datasets.visium). The notebook will preprocess the spatial data and generate either single map or double map. If you wish to compare two single cell data, you should filter at **3.3.1 Separate to two single cell data**.

### Output
If you are using **Tangram_CPU**, you should have 2 output files: 
- adata_st
- ad_map

If you are using **Tangram_CPU_Double**, you should have 3 output files: 
- adata_st
- ad_map1
- ad_map2

## Tangram_CPU
### Input
**Tangram_CPU** needs 3 input files: 
- adata_sc: pasca.log1p_liger_med_singleR_noglyc.h5ad
- adata_st
- ad_map

**Tangram_CPU_Double** needs 4 input files:
- adata_sc: pasca.log1p_liger_med_singleR_noglyc.h5ad
- adata_st
- ad_map1
- ad_map2

### Checklist

1. **Filter values** must be the same for both Tangram_GPU and Tangram_CPU
2. **Marker genes** must be the same for both Tangram_GPU and Tangram_CPU

### Process
The choice of genes to visualize can be done in the variable **GOI** in **2.4 - inspect predictions**.


### Output
Tangram_CPU pipeline has 0 output files.
