import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy as sp
from scipy.sparse import csr_matrix
import pandas as pd
import seaborn as sns
import pathlib
import git

import rpy2.rinterface_lib.callbacks
import rpy2.robjects as ro
import logging
from rpy2.robjects import numpy2ri
from rpy2.robjects import pandas2ri
import anndata2ri

pandas2ri.activate()
anndata2ri.activate()

git_repo = git.Repo(".", search_parent_directories=True)
git_root = git_repo.git.rev_parse("--show-toplevel")
topDir = os.path.join(git_root, "raw_data/Pasca_scRNAseq/")
prefix = 'pasca'
pathlib.Path(topDir)

inDir = pathlib.Path(topDir,"raw")
outDir = pathlib.Path(topDir,"outputs")
interDir = pathlib.Path(topDir,"inter")

def fullpath_closure(data_dir: pathlib.Path) -> pathlib.Path | str:
    import os
    def fullpath(path):
        return os.path.join(data_dir, path)
    return fullpath

def concatenate(directories: list[str]) -> sc.AnnData:
    """
    Concatenates scanpy annotation data from multiple samples into 1 large object
    """
    # Each directory corresponds with a sample
    raw_data = [sc.read(dir) for dir in directories]
    
    # Take the first data object, and concatenate with the remaining objects
    concatenated = raw_data[0].concatenate(*raw_data[1:])
    return concatenated

def sparsify(anndata: sc.AnnData) -> sc.AnnData:
    # Ensures we're not storing redundant variables in data
    anndata.var_names_make_unique()

    # CSR is compressed sparse row format; simply a space saving measure
    anndata.X = csr_matrix(anndata.X)
    return anndata 


def annotate(anndata: sc.AnnData, directories: list[str]) -> None:
    """
    Annotates scanpy AnnData object with condition and sample IDs from file names
    """

    # Get file names in case full paths are passed in
    basenames = [os.path.basename(d) for d in directories]

    # Check that data is labeled in a way expected by remaining logic
    try:
        for dir in basenames:
            cond = dir.split("_")[-2].lower()
            assert (cond == 'control') or (cond == 'patient')
    except AssertionError:
        print(f"This path name doesn't contain the correct condition label in the correct location: {dir}")
        print(f"Please rename the path")
        sys.exit(255)

    # Find each experimental condition (usually control and patient)
    condition_map = {i: basenames[i].split("_")[-2] for i in range(len(basenames))}
    sample_id_map = {i: basenames for i in range(len(basenames))}

    # Map corresponding variables to respective dictionaries
    anndata.obs['condition'] = anndata.obs['batch'].map(condition_map)
    anndata.obs['sampleID'] = anndata.obs['batch'].map(sample_id_map)

def write_anndata(anndata: sc.AnnData, dir: pathlib.Path, filename: str) -> None:
    """
    Writes scanpy AnnData object to disk in specified directory with filename
    """
    anndata.write(pathlib.Path(dir, filename))

def adata_var_get(_adata: sc.AnnData, prefix_l=[], gene_l=[]):
    vars_l = set()
    for prefix in prefix_l:
        vars_l |= set(_adata.var.index[_adata.var.index.str.startswith(prefix)])
    vars_l |= set(_adata.var.index[_adata.var.index.isin(gene_l)])
    if vars_l:
        return _adata[:, list(vars_l)]
    else:
        return None

def highest_expr_genes(anndata: sc.AnnData, **kwargs) -> plt.Figure:
    """
    Plots highest expressed genes from scanpy AnnData object
    """
    
    return sc.pl.highest_expr_genes(anndata, **kwargs)

def mito_qc_metrics(anndata: sc.AnnData) -> pd.DataFrame:
    """
    Uses mitochondrial genes to calculate quality control metrics for annotation data
    """

    anndata.var['mt'] = anndata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    qc_metrics = sc.pp.calculate_qc_metrics(
        anndata, 
        qc_vars=['mt'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    return qc_metrics



if __name__ == "__main__":
    # Within the raw folder, get all the directory paths
    anndata_dirs = [p for p in os.listdir(inDir) if os.path.isdir(p)]
    
    # Concatenate all these files into one scanpy object
    anndata_concat = concatenate(anndata_dirs)

    # Save space by storing sparse matrices in CSR format
    anndata_sparse = sparsify(anndata_concat)

    # Store the raw object in case
    write_anndata(anndata_sparse, interDir, "pasca_aggr.h5ad")

    # Annotate control/patient labels according to directory names
    anndata_annot = annotate(anndata_sparse, anndata_dirs)
    
    ribo_p = adata_var_get(anndata_annot, prefix_l=['RPL', 'RPS', 'MRPL', 'MRPS'])
    ribo_r = adata_var_get(anndata_annot, gene_l=['RN45S', 'RN4.5S'])
    print(ribo_p)
    print(ribo_r)