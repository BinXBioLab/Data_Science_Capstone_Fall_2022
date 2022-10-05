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
    Annotates scanpy AnnData objects with appropriate condition and sample IDs
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



