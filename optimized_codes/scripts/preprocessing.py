import pdb
import gc
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import pandas as pd
import pathlib
import git
from tqdm import tqdm, trange
import click

import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import anndata2ri

plt.ioff()            # Ensures non-interactive matplotlib to avoid pop-up figure windows 
pandas2ri.activate()  # Pandas interface from python <--> R
anndata2ri.activate() # AnnData interface from python <--> R

# TODO: Make this so it recognizes Google Colab versus local machine
# TODO: Add logic to switch between looking for GitHub root and colab
git_repo = git.Repo(".", search_parent_directories=True)
git_root = git_repo.git.rev_parse("--show-toplevel")

topDir = os.path.join(git_root, "raw_data/Pasca_scRNAseq/")
prefix = 'pasca'
inDir = pathlib.Path(topDir,"raw")
outDir = pathlib.Path(topDir,"outputs")
interDir = pathlib.Path(topDir,"inter-test")

def fullpath_closure(data_dir: pathlib.Path) -> pathlib.Path | str:
    import os
    def fullpath(path):
        return os.path.join(data_dir, path)
    return fullpath

def concatenate(directories: list[str]) -> sc.AnnData:
    """
    Concatenates scanpy annotation data from multiple samples into 1 large object
    """
    click.echo(f"Reading in 10x Genomics counts and metadata")

    # Each directory corresponds with a sample
    raw_data = [sc.read_10x_mtx(dir, var_names='gene_symbols', cache=True, gex_only=True) for dir in tqdm(directories)]
    
    # Take the first data object, and concatenate with the remaining objects
    concatenated = raw_data[0].concatenate(*raw_data[1:])
    return concatenated

def sparsify(anndata: sc.AnnData, make_vars_unique: bool = True) -> sc.AnnData:
    click.echo(f"Sparsifying AnnData counts matrix")
    
    # Ensures we're not storing redundant variables in data
    if make_vars_unique:
        anndata.var_names_make_unique()

    # CSR is compressed sparse row format; simply a space saving measure
    anndata.X = csr_matrix(anndata.X)
    return anndata 

def annotate(anndata: sc.AnnData, directories: list[str]) -> None:
    """
    Annotates scanpy AnnData object with condition and sample IDs from file names
    """

    click.echo(f"Annotating AnnData object")
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
    return anndata

def write_anndata(anndata: sc.AnnData, dir: pathlib.Path, filename: str) -> None:
    """
    Writes scanpy AnnData object to disk in specified directory with filename
    """
    fullpath = pathlib.Path(dir, filename)
    click.echo(f"Writing AnnData object to {fullpath}")
    anndata.write(fullpath)

def adata_var_get(anndata: sc.AnnData, prefix_l=[], gene_l=[]):
    vars_l = set()
    for prefix in prefix_l:
        vars_l |= set(anndata.var.index[anndata.var.index.str.startswith(prefix)])
    vars_l |= set(anndata.var.index[anndata.var.index.isin(gene_l)])
    if vars_l:
        return anndata[:, list(vars_l)]
    else:
        return None

def mito_qc_metrics(anndata: sc.AnnData) -> pd.DataFrame:
    """
    Uses mitochondrial genes to calculate quality control metrics for annotation data
    """

    anndata.var['mt'] = anndata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(anndata, qc_vars=['mt'], percent_top=None, log1p=False,  inplace=True)
    return anndata

def augment_with_ribo_data(anndata: sc.AnnData, ribo_p: any) -> None:
    anndata.obs['n_counts_ribo_p'] = ribo_p.X.sum(axis=1).A1.tolist()

    # TODO: Turn this into multiple functions for improved readability
    anndata.obs['percent_ribo_p'] = (100*ribo_p.X.sum(axis=1).transpose()/np.ravel(np.sum(anndata.X, axis=1))).transpose().A1.tolist()

    return anndata

def create_qc_violin_plots(anndata: sc.AnnData, to_plot):
    progress = trange(3, desc='Creating quality control violin plots')

    # Basic QC
    progress.set_description("Basic QC plot")
    sc.pl.violin(anndata, to_plot, jitter=0.4, multi_panel=True, save='.qc.pdf')
    progress.update(1)

    # Grouped by batch
    progress.set_description("Grouped by batch")
    sc.pl.violin(anndata, to_plot, jitter=0.4, save=f'.qc.batch.pdf', groupby='batch', rotation=45, cut=0)
    progress.update(1)
    
    # Grouped by condition
    progress.set_description("Grouped by condition")
    sc.pl.violin(anndata, to_plot, jitter=0.4, save=f'.qc.batch.pdf', groupby='condition', rotation=45, cut=0)
    progress.update(1)

def filter_by_attr(anndata: sc.AnnData,  pct_count, genes_by_count, total_counts) -> sc.AnnData:
    # Perform a deep copy so operations are not done in-place
    temp = anndata.copy()

    # Filter by percentage of MT- genes
    if pct_count[0] is not None:
        temp = temp[temp.obs['pct_counts_mt'] > pct_count[0], :]
    if pct_count[1] is not None:
        temp = temp[temp.obs['pct_counts_mt'] < pct_count[1], :]
    
    # Filter by number of detected genes
    if genes_by_count[0] is not None:
        temp = temp[temp.obs['n_genes_by_counts'] >= genes_by_count[0], :]
    if genes_by_count[1] is not None:
        temp = temp[temp.obs['n_genes_by_counts'] <= genes_by_count[1], :]
    
    # Filter by number of counts
    if total_counts[0] is not None:
        temp = temp[temp.obs['total_counts'] >= total_counts[0], :]
    if total_counts[1] is not None:
        temp = temp[temp.obs['total_counts'] <= total_counts[1], :]

    return temp

def r_pipeline(anndata: sc.AnnData):
    # Importing R package scran
    importr('scran')
    
    # Assigning the transposed counts to a variable called `mat`
    ro.r.assign('mat', anndata.X.T)
    qclust_params = 'mat'

    # Can this be done in python?
    # Creating quickCluster object and transforming into SingleCellExperiment for later use
    ro.reval(f'cl <- quickCluster({qclust_params})')
    csf_params = f'SingleCellExperiment::SingleCellExperiment(list(counts=mat)), clusters = cl'

    # Garbage collection to free up memory
    del qclust_params
    gc.collect()

    # Something in scran?
    ro.reval(f'mysce <-computeSumFactors({csf_params})')

    # Converting scran output to numpy
    sf = np.asarray(ro.reval(f'sizeFactors(mysce)'))

    # Garbage collection again?
    del csf_params
    gc.collect()

    # TODO: Change this so manipulation doesn't occur in-place
    anndata.obs['sf'] = sf   
    anndata.layers['counts'] = anndata.X.copy()
    anndata.X /= anndata.obs['sf'].values[:, None]

    return anndata

def generate_liger_input(anndata: sc.AnnData) -> None:
    # Counts file needed by liger
    anndata.X = anndata.layers['counts']
    # anndata.to_df().to_csv(pathlib.Path(interDir, 'pasca_log1p.csv'))
    anndata.to_df().to_parquet(pathlib.Path(interDir, 'pasca_log1p.csv')) # Using parquet instead of csv to save space

    # Metadata file needed by liger
    sc.get.obs_df(anndata, keys=anndata.obs_keys()).to_csv(
        pathlib.Path(interDir,'pasca_log1p.metadata.csv')
    )

def preprocess():
    # Within the raw folder, get all the directory paths
    anndata_dirs = [os.path.join(inDir, p) for p in os.listdir(inDir) if os.path.isdir(os.path.join(inDir, p))]
    # pdb.set_trace()
    
    # Concatenate all these files into one scanpy object
    anndata_concat = concatenate(anndata_dirs)

    # Save space by storing sparse matrices in CSR format
    anndata_sparse = sparsify(anndata_concat, make_vars_unique=True)

    # Store the raw object in case
    write_anndata(anndata_sparse, interDir, "pasca_aggr.h5ad")

    # Annotate control/patient labels according to directory names
    anndata_annot = annotate(anndata_sparse, anndata_dirs)
    
    # Create violin plot of highest expressed genes
    sc.pl.highest_expr_genes(anndata_annot, n_top=40, save='.raw_top40.pdf');
    
    # TODO: How to not have this hardcoded
    # Get ribosomal gene data
    ribo_p = adata_var_get(anndata_annot, prefix_l=['RPL', 'RPS', 'MRPL', 'MRPS'])
    ribo_r = adata_var_get(anndata_annot, gene_l=['RN45S', 'RN4.5S'])
    print(f"ribo_p: {ribo_p}")
    print(f"ribo_r: {ribo_r}")

    anndata_augmented = mito_qc_metrics(anndata_annot)
    anndata_augmented = augment_with_ribo_data(anndata_annot, ribo_p)
    items_to_plot = [
        'n_genes_by_counts',
        'total_counts',
        'pct_counts_mt',
        'n_counts_ribo_p',
        'percent_ribo_p',
    ]

    # Creating violin plots for visual quality control
    create_qc_violin_plots(anndata=anndata_augmented, to_plot=items_to_plot)

    # Save anndata before QC step
    write_anndata(anndata_annot, interDir, "pasca_preqc.h5ad")

    # Plots before filtering
    sc.pl.scatter(
        anndata_annot, 
        color='pct_counts_mt', 
        x='total_counts', 
        y='n_genes_by_counts', 
        save='.qc.lib_detect_mito.pdf', 
        color_map='viridis'
    )
    sc.pl.scatter(
        anndata_annot, 
        color='percent_ribo_p', 
        x='total_counts', 
        y='n_genes_by_counts', 
        save='.qc.lib_detect_ribo_p.pdf', 
        color_map='viridis'
    )
    
    # Filtering
    anndata_filtered = filter_by_attr(anndata=anndata_annot, pct_count=(None, 20), genes_by_count=(200, 8000), total_counts=(None, 5e4))

    # Plots after filtering
    sc.pl.scatter(
        anndata_filtered, 
        color='pct_counts_mt', 
        x='total_counts', 
        y='n_genes_by_counts', 
        frameon=True, 
        save='.qc.lib_detect_mito.filter.pdf', 
        color_map='viridis'
    )
    sc.pl.scatter(
        anndata_filtered,
        color='percent_ribo_p',
        x='total_counts',
        y='n_genes_by_counts',
        frameon=True,
        save='.qc.lib_detect_ribo_p.filter.pdf',
        color_map='viridis',
    )
    sc.pl.scatter(
        anndata_filtered,
        color='batch',
        x='total_counts',
        y='n_genes_by_counts',
        frameon=True,
        save='.qc.lib_detect_batch.filter.pdf',
        color_map='viridis',
    )

    # Filter out genes with 0 cells
    sc.pp.filter_genes(anndata_filtered, min_cells=1)

    # Save anndata before QC step
    write_anndata(anndata_filtered, interDir, "pasca_postqc.h5ad")

    # scran pipeline in R
    anndata_scran = r_pipeline(anndata_filtered)

    # Sparsifying the scran output to save space
    anndata_scran_sparse = sparsify(anndata=anndata_scran, make_vars_unique=False)

    # Plotting cell counts
    sc.pl.scatter(
        anndata_scran_sparse,
        color='total_counts',
        y='sf',
        x='n_genes_by_counts',
        color_map='viridis',
        save='.sf.pdf',
    )

    # Take each entry in the counts matrix (call it x) and replace it with log(x + 1)
    sc.pp.log1p(anndata_scran_sparse)
    write_anndata(anndata_scran_sparse, interDir, "pasca_log1p.h5ad")

if __name__ == "__main__":
    preprocess()
