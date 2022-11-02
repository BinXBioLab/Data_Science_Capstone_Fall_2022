import gc
import scanpy as sc
import git
import pathlib
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from helper import write_anndata
import wget

git_repo = git.Repo(".", search_parent_directories=True)
git_root = git_repo.git.rev_parse("--show-toplevel")

topDir = os.path.join(git_root, "raw_data/Pasca_scRNAseq/")
prefix = 'pasca'
inDir = pathlib.Path(topDir,"raw")
outDir = pathlib.Path(topDir,"outputs")
interDir = pathlib.Path(topDir,"inter-test")

sns.set_style('white')
plt.rcParams['savefig.facecolor'] = 'w'
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.set_figure_params(scanpy=True, dpi=120, dpi_save=150, frameon=True, vector_friendly=True, fontsize=14) 
sc.settings.figdir = outDir('fig_supp')

# TODO: Check later if this can be removed
samples = ['pasca']
prefix = samples[0]
prefix

def generate_figures(anndata: sc.AnnData) -> None:
    # Compute neighborhood graph
    # use_rep means use X_liger result to represent X for calculation
    # other adjustable parameters that control uMAP output are n_neighbors, metrics(distance), min_dist
    sc.pp.neighbors(
        anndata,
        n_neighbors=15,
        random_state=42,
        use_rep= 'X_liger',
    )

    # conduct umap
    sc.tl.umap(
        anndata, 
        random_state=42, 
        spread=3.0, 
        min_dist=0.1
    )

    # Plot UMAP
    sc.pl.umap(anndata, color=["PAX6"], show=False) 
    sc.tl.leiden(anndata)
    sc.pl.umap(anndata, color=['leiden'], show=False)

    return anndata

def rank_genes_and_save(anndata: sc.AnnData) -> None:
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    write_anndata(anndata, interDir, "pasca_log1p_liger.h5ad")

def download_and_process_reference(
    counts_url=None, 
    counts_name=None, 
    meta_url=None, 
    meta_name=None) -> sc.AnnData:
    
    # Download reference material
    if counts_url == None:
        matrix_filepath = os.path.join(interDir, "exprMatrix.tsv.gz")
        wget.download("https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz", matrix_filepath)
    else:
        matrix_filepath = os.path.join(interDir, counts_name)
        wget.download(counts_url, matrix_filepath)

    if meta_url == None:
        meta_filepath = os.path.join(interDir, "meta.tsv")
        wget.download("https://cells.ucsc.edu/cortex-dev/meta.tsv", meta_filepath)
    else:
        meta_filepath = os.path.join(interDir, meta_name)
        wget.download(meta_url, meta_filepath)

    # Combine counts and metadata
    h_ad = sc.read_text(matrix_filepath).transpose()
    meta = pd.read_csv(meta_filepath)
    h_ad.obs = meta
    anndata = sc.AnnData(h_ad)

    return anndata



if __name__ == "__main__":
    # Read files
    path_to_anndata = os.path.join(interDir, "pasca_log1p.h5ad")
    path_to_liger = os.path.join(interDir, "pasca.liger.csv")
    df_liger = pd.read_csv(path_to_liger, index_col=0)
    adata = sc.read(path_to_anndata)

    # Generating figures based on Liger output
    adata.obsm['X_liger'] = df_liger.loc[adata.obs_names, :].values
    adata = generate_figures(adata)
    rank_genes_and_save(adata)

    # ???
    postqc_path = os.path.join(interDir, "pasca_postqc.h5ad")
    adata = sc.read(postqc_path)
    df = pd.DataFrame.sparse.from_spmatrix(
        adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )

    # filter genes?
    df = df.loc[:, ((df > 0).sum(axis=0) >= 1)]

    ad_nowakowski = download_and_process_reference()
    sc.pp.filter_cells(ad_nowakowski, min_genes=200)
    sc.pp.filter_genes(ad_nowakowski, min_cells=1)
    genes = df.columns.intersection(ad_nowakowski.var_names)

    df = df.loc[:, genes]
    ad_nowakowski = ad_nowakowski[:, genes].copy()

    # df is the target dataframe to be exported to R 
    # generate the ref dataframe with CellID as index to be exported to R
    df_nowakowski = pd.DataFrame(
        ad_nowakowski.X,
        index=ad_nowakowski.obs['Cell'],
        columns=ad_nowakowski.var_names
    )

    # get the annotation (WGCNAcluster in this case) for each CellID (Cell in this case)
    ad_nowakowski_label = pd.DataFrame(
        ad_nowakowski.obs["WGCNAcluster"].values,
        index=ad_nowakowski.obs['Cell']
    )
    ad_nowakowski_label.columns = ['WGCNAcluster']
    ad_nowakowski_label.to_csv(os.path.join(interDir, 'nowakowski.label_bx.csv'), sep=',', header=True)
    
    # Creating noglyc version of nowakowski?
    ad_nowakowski_noglyc = ad_nowakowski[
        ~(ad_nowakowski.obs['WGCNAcluster'] == 'Glyc'),
        :
    ].copy()
    ad_nowakowski_noglyc_label = pd.DataFrame(ad_nowakowski_noglyc.obs)
    ad_nowakowski_noglyc_label_f = ad_nowakowski_noglyc_label[["WGCNAcluster"]]
    ad_nowakowski_noglyc_label_f.index = ad_nowakowski_noglyc_label["Cell"]
    ad_nowakowski_noglyc_label_f.to_csv(
        os.path.join(interDir, 
        'nowakowski_noglyc.label_bx.csv'), 
        sep=',', 
        header=True
    )






