import pdb
from gc import collect
import sys
import os
import pathlib
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from git import Repo
from git.exc import InvalidGitRepositoryError
from tqdm import tqdm, trange
import click

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import anndata2ri

sc.settings.verbosity = 0  # Only print ERROR level messages to the console
pandas2ri.activate()       # Pandas interface from python <--> R
anndata2ri.activate()      # AnnData interface from python <--> R
prefix = "pasca"


def write_anndata(anndata: sc.AnnData, directory: pathlib.Path,
                  filename: str) -> None:
    """
    Writes scanpy AnnData object to disk in specified directory with filename
    """
    fullpath = pathlib.Path(directory, filename)
    click.echo(f"\nWriting AnnData object to {fullpath}\n")
    anndata.write(fullpath)


def concatenate(directories: list, cache=False) -> sc.AnnData:
    """
    Concatenates scanpy annotation data from multiple samples into 1 large object

    Args:
        directories (list): List of folder locations for 10x Genomic data

    Returns:
        sc.AnnData: All controls and patients concatenated
    """
    click.echo("Reading in 10x Genomics counts and metadata")

    # Each directory corresponds with a sample
    # WARNING: AnnData performs implicit np.float64 -> np.float32 conversion, but this will break in 0.9
    # Sorted order is important to do all controls followed by all patients
    # (might need to be refactored in the future)
    loop = [item for item in directories if not item.endswith('cache')] # ignore cache directory
    loop = tqdm(sorted(loop))
    raw_data = []

    for folder in loop:
        loop.set_postfix(directory=os.path.basename(folder))
        sample = sc.read_10x_mtx(
            folder,
            cache=cache,
            var_names='gene_symbols',
            gex_only=True
        )
        raw_data.append(sample)

    # Take the first data object, and concatenate with the remaining objects
    concatenated = raw_data[0].concatenate(*raw_data[1:])
    return concatenated


def sparsify(anndata: sc.AnnData, make_vars_unique: bool = True) -> sc.AnnData:
    """
    Transform counts matrix into compressed sparse row (CSR) format to save memory/disk space

    Args:
        anndata (sc.AnnData): Annotated data object to sparsify
        make_vars_unique (bool): Whether or not to make the variable names unique (default True)

    Returns:
        sc.AnnData: Annotated data with sparsified counts
    """
    click.echo("\nSparsifying AnnData counts matrix\n")

    # Ensures we're not storing redundant variables in data
    if make_vars_unique:
        anndata.var_names_make_unique()

    # CSR is compressed sparse row format; simply a space saving measure
    anndata.X = csr_matrix(anndata.X)
    return anndata


def _annotate(anndata: sc.AnnData, directories: list) -> None:
    """
    Internal function to annotate object with condition and sample IDs from file names

    Args:
        anndata (sc.AnnData): Annotated data object
        directories (list): List of 10x Genomics folders (each with matrix, metadata, labels)
    """

    click.echo("Annotating AnnData object")
    # Get file names in case full paths are passed in
    basenames = [os.path.basename(d) for d in directories]

    # Check that data is labeled in a way expected by remaining logic
    try:
        for directory in basenames:
            condition = directory.split("_")[-2].lower()
            assert (condition == 'control') or (condition == 'patient')
    except AssertionError:
        print(f"This path name doesn't contain the correct condition label in the correct location: {directory}")
        print("Please rename the path")
        sys.exit(255)

    # Find each experimental condition (usually control and patient)
    condition_map = {str(i): basenames[i].split("_")[-2]
                     for i in range(len(basenames))}
    sample_id_map = {
        str(idx): basename for idx,
        basename in enumerate(basenames)}

    # Map corresponding variables to respective dictionaries
    anndata.obs['condition'] = anndata.obs['batch'].map(condition_map)
    anndata.obs['sampleID'] = anndata.obs['batch'].map(sample_id_map)
    return anndata


def adata_var_get(anndata: sc.AnnData, prefix_l=[], gene_l=[]):
    """
    Get unique variables?

    Args:
        anndata (sc.AnnData): Annotated data object
        prefix_l (list): ???
        gene_l (list): ???

    Returns:
        sc.AnnData: Annotated data object
    """
    vars_l = set()

    for prefix in prefix_l:
        vars_l |= set(
            anndata.var.index[anndata.var.index.str.startswith(prefix)])
    vars_l |= set(anndata.var.index[anndata.var.index.isin(gene_l)])

    return anndata[:, list(vars_l)] if vars_l else None


def mito_qc_metrics(anndata: sc.AnnData) -> sc.AnnData:
    """
    Uses mitochondrial genes to calculate quality control metrics for annotation data

    Args:
        anndata (sc.AnnData): Annotated data object

    Returns:
        sc.AnnData: Same annotated data object as input
    """

    # annotate the group of mitochondrial genes as 'mt'
    anndata.var['mt'] = anndata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        anndata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    return anndata


def augment_with_ribo_data(anndata: sc.AnnData, ribo_p: any) -> None:
    """
    Add ribosomal proteins

    Args:
        anndata (sc.AnnData): Annotated data object

    Returns:
        sc.AnnData: Same annotated data object as input
    """
    anndata.obs['n_counts_ribo_p'] = ribo_p.X.sum(axis=1).A1.tolist()

    # TODO: Turn this into multiple functions for improved readability
    anndata.obs['percent_ribo_p'] = (
        100 *
        ribo_p.X.sum(
            axis=1).transpose() /
        np.ravel(
            np.sum(
                anndata.X,
                axis=1))).transpose().A1.tolist()

    return anndata


def create_qc_violin_plots(anndata: sc.AnnData, to_plot: list) -> None:
    """
    Create violin plots for the purposes of visual quality control

    Args:
        anndata (sc.AnnData): Annotated data object
        to_plot (list): List of genes to plot
    """

    progress = trange(3, desc='Creating quality control violin plots')

    # Basic QC
    progress.set_description("Basic QC plot")
    sc.pl.violin(
        anndata,
        to_plot,
        jitter=0.4,
        multi_panel=True,
        save='.qc.pdf',
        show=False)
    progress.update(1)

    # Grouped by batch
    progress.set_description("Grouped by batch")
    sc.pl.violin(
        anndata,
        to_plot,
        jitter=0.4,
        save=f'.qc.batch.pdf',
        groupby='batch',
        rotation=45,
        cut=0,
        show=False)
    progress.update(1)

    # Grouped by condition
    progress.set_description("Grouped by condition")
    sc.pl.violin(
        anndata,
        to_plot,
        jitter=0.4,
        save=f'.qc.batch.pdf',
        groupby='condition',
        rotation=45,
        cut=0,
        show=False)
    progress.update(1)


def filter_by_attr(anndata: sc.AnnData, pct_count,
                   genes_by_count, total_counts) -> sc.AnnData:
    """
    Filter counts matrix based on percent counts, number of genes, and total cell counts

    Args:
        anndata (sc.AnnData): Annotated data object
        pct_count (tuple): Bounds for percent count filtering of the form (lower bound, upper bound)
        genes_by_count (tuple): Bounds for filtering low or high gene counts of the form (lower bound, upper bound)
        total_counts (tuple): Bounds for filtering low and/or high total counts of the form (lower bound, upper bound)

    Returns:
        sc.AnnData: Annotated data object
    """
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


def r_pipeline(anndata: sc.AnnData, inter_dir: str):
    """
    Run scran and quickClister on anndata using R

    Args:
        anndata (sc.AnnData): Annotated data object
    """
    try:
        click.echo("Running SCRAN as Script\n")
        from scipy.io import mmwrite
        
        input_path = os.path.join(inter_dir, "scran_input.mtx")
        output_path = os.path.join(inter_dir, "scran_output.mtx")
        print(f"Input path: {input_path}")
        print(f"Output path: {output_path}")
        
        mmwrite(input_path, anndata.X.T)
        cmd = f"Rscript --vanilla ../R/run_scran.R {input_path} {output_path}"
        print(cmd)
        os.system(cmd)
        
        sf = np.genfromtxt(output_path)
        anndata.obs['sf'] = sf
        anndata.layers['counts'] = anndata.X.copy()
        anndata.X /= anndata.obs['sf'].values[:, None]
    except FileNotFoundError:
        click.echo("Couldn't find one or more of the files for scran")
        click.echo("Falling back to rpy2")
        progress = trange(6, desc='R single cell experiment pipeline')

        # Importing R package scran
        progress.set_description("Importing scran package")
        importr('scran')
        progress.update(1)

        # Assigning the transposed counts to a variable called `mat`
        progress.set_description("Storing anndata matrix in R object")
        ro.r.assign('mat', anndata.X.T)
        progress.update(1)

        # TODO: Can this be done in python?
        # Creating quickCluster object and transforming into SingleCellExperiment
        # for later use
        progress.set_description("quickCluster")
        ro.reval('cl <- quickCluster(mat)')
        csf_params = 'SingleCellExperiment::SingleCellExperiment(list(counts=mat)), clusters = cl'
        progress.update(1)

        # Garbage collection to free up memory
        # del qclust_params
        collect()

        # Something in scran?
        progress.set_description("computeSumFactors")
        ro.reval(f'mysce <-computeSumFactors({csf_params})')
        progress.update(1)

        # Converting scran output to numpy
        progress.set_description("sizeFactors")
        sf = np.asarray(ro.reval('sizeFactors(mysce)'))
        progress.update(1)

        # Garbage collection again?
        del csf_params
        collect()

        # TODO: Change this so manipulation doesn't occur in-place
        progress.set_description("Assigning R results to python")
        anndata.obs['sf'] = sf
        anndata.layers['counts'] = anndata.X.copy()
        anndata.X /= anndata.obs['sf'].values[:, None]
        progress.update(1)

    return anndata


def generate_liger_input(anndata: sc.AnnData, inter_dir) -> None:
    """
    Create appropriate counts and metadata files that can be used in LigeR

    Args:
        anndata (sc.AnnData): Annotated data object
    """
    # Counts file needed by liger
    anndata.X = anndata.layers['counts']

    # Using parquet instead of csv to save space
    anndata.to_df().to_parquet(pathlib.Path(inter_dir, 'pasca_log1p.parquet'))

    # Metadata file needed by liger
    sc.get.obs_df(anndata, keys=anndata.obs_keys()).to_csv(
        pathlib.Path(inter_dir, 'pasca_log1p.metadata.csv'))


def aggregate(root, aggregated_filename, inter_dir):
    """High level function which combines loading, concatenating, and sparsifying the data

    Args:
        root (str): Root directory where 10x Genomics files are stored
        aggregated_filename (str): Name to give the aggregated anndata object when saving to disk

    Returns:
        list: List of directories containing original counts data files
        sc.AnnData: Sparsified, annotated data
    """
    # Within the raw folder, get all directory paths as absolute paths
    anndata_dirs = [
        os.path.join(root, p)
        for p in os.listdir(root)
        if os.path.isdir(os.path.join(root, p))
    ]

    # Concatenate all these counts and labels into one AnnData object
    anndata_concat = concatenate(anndata_dirs)

    # Save space by storing sparse matrices in CSR format
    anndata_sparse = sparsify(anndata_concat, make_vars_unique=True)

    # Store the unprocessed, concatenated data
    write_anndata(anndata_sparse, inter_dir, aggregated_filename)
    return anndata_dirs, anndata_sparse


def annotate(anndata_sparse: sc.AnnData, anndata_dirs: list, inter_dir):
    """
    High-level function which combines annotation, getting unique variables, and augmenting with ribosomal data

    Args:
        anndata_sparse (sc.AnnData): Sparsified annotated data object
        anndata_dirs (list): List of directories containing original counts data files

    Returns:
        sc.AnnData: Returns annodated data that's been augmented with info about ribosomal proteins
    """
    # Annotate control/patient labels according to directory names
    anndata_annot = _annotate(anndata_sparse, anndata_dirs)

    # Create violin plot of highest expressed genes
    sc.pl.highest_expr_genes(
        anndata_annot,
        n_top=40,
        save='.raw_top40.pdf',
        show=False)

    # TODO: How to not have this hardcoded
    # Get ribosomal gene data
    ribo_p = adata_var_get(anndata_annot, prefix_l=['RPL', 'RPS', 'MRPL', 'MRPS'])

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
    write_anndata(anndata_annot, inter_dir, "pasca_preqc.h5ad")
    return anndata_annot


def quality_control(anndata_annot: sc.AnnData, inter_dir, **kwargs):
    """High-level function which performs all quality control steps including plotting and filtering

    Args:
        anndata_annot (sc.AnnDat): Annotated data that's been augmented with ribosomal proteins
        **kwargs (dict): Additional keyword arguments to pass to attribute filtering function

    Returns:
        sc.AnnData: Annotated data object after filtering by percent counts, gene counts, and total counts
    """
    # Plots before filtering
    scatter_kwargs = {
        'x': 'total_counts',
        'y': 'n_genes_by_counts',
        'color_map': 'viridis',
        'show': False
    }
    sc.pl.scatter(
        anndata_annot,
        color='pct_counts_mt',
        save='.qc.lib_detect_mito.pdf',
        **scatter_kwargs)
    sc.pl.scatter(
        anndata_annot,
        color='percent_ribo_p',
        save='.qc.lib_detect_ribo_p.pdf',
        **scatter_kwargs)

    # Filtering
    anndata_filtered = filter_by_attr(anndata=anndata_annot, **kwargs)

    # Plots after filtering
    sc.pl.scatter(
        anndata_filtered,
        color='pct_counts_mt',
        frameon=True,
        save='.qc.lib_detect_mito.filter.pdf',
        **scatter_kwargs)
    sc.pl.scatter(
        anndata_filtered,
        color='percent_ribo_p',
        frameon=True,
        save='.qc.lib_detect_ribo_p.filter.pdf',
        **scatter_kwargs)
    sc.pl.scatter(
        anndata_filtered,
        color='batch',
        frameon=True,
        save='.qc.lib_detect_batch.filter.pdf',
        **scatter_kwargs)

    # Filter out genes with 0 cells
    sc.pp.filter_genes(anndata_filtered, min_cells=1, inplace=True, copy=True)

    # Save anndata before QC step
    write_anndata(anndata_filtered, inter_dir, "pasca_postqc.h5ad")
    return anndata_filtered


def postprocess(anndata_filtered: sc.AnnData, inter_dir):
    """
    Final steps in preprocessing pipeline including running scran, quickCluster, and log transforming counts

    Args:
        anndata_filtered (sc.AnnData): Annotated data object after filtering by percent counts, gene counts, and total counts
    """
    # scran pipeline in R
    anndata_scran = r_pipeline(anndata_filtered, inter_dir)

    # Sparsifying the scran output to save space
    anndata_scran_sparse = sparsify(
        anndata=anndata_scran,
        make_vars_unique=False)

    # Plotting cell counts
    sc.pl.scatter(
        anndata_scran_sparse,
        color='total_counts',
        y='sf',
        x='n_genes_by_counts',
        color_map='viridis',
        save='.sf.pdf',
        show=False
    )

    # Take each entry in the counts matrix (call it x) and replace it with log(x + 1)
    sc.pp.log1p(anndata_scran_sparse)
    write_anndata(anndata_scran_sparse, inter_dir, "pasca_log1p.h5ad")
    generate_liger_input(anndata_scran_sparse, inter_dir)


@click.command()
@click.option('--in_dir', default=None, help="Parent directory of 10x Genomics directories")
@click.option('--inter_dir', default=None, help="Directory to store intermediate outputs")
@click.option('--pct_count_lower', default=None, help='Lower bound of percent counts to filter')
@click.option('--pct_count_upper', default=20, help='Upper bound of percent counts to filter')
@click.option('--genes_by_count_lower', default=200, help='Lower bound of genes by counts to filter')
@click.option('--genes_by_count_upper', default=8000, help='Upper bound of genes by counts to filter')
@click.option('--total_counts_lower', default=None, help='Lower bound of total counts to filter')
@click.option('--total_counts_upper', default=5e4, help='Upper bound of total counts to filter')
@click.option('--aggregated_filename', default="pasca_aggr.h5ad", help='Filename for concatenated counts + labels h5ad data file')
def preprocess(
        in_dir: str,
        inter_dir: str,
        pct_count_lower: int,
        pct_count_upper: int,
        genes_by_count_lower: int,
        genes_by_count_upper: int,
        total_counts_lower: float,
        total_counts_upper: float,
        aggregated_filename) -> None:
    """Top-level function which executes all the high-level functions"""
    
    # Try automatically finding the repo root if the inputs are not specified
    if in_dir is None:
        try:
            git_repo = Repo(".", search_parent_directories=True)
            git_root = git_repo.git.rev_parse("--show-toplevel")

            in_dir = pathlib.Path(git_root, "raw_data/Pasca_scRNAseq/", "raw")
            inter_dir = pathlib.Path(git_root, "raw_data/Pasca_scRNAseq/", "inter-test")
        except InvalidGitRepositoryError:
            error = """
            Looks like there's no git repository :(

            To avoid this error, run the preprocessing script with 
                --in_dir set to the directory containing your inputs
                --inter_dir set to the directory you want to store the intermediate files

            Example: python preprocessing.py --in_dir ./path/to/inputs --inter_dir ./inter/ ...
            """
            print(error)
    else:
        in_dir = pathlib.Path(in_dir)
    
    # Store the intermediate files in the current working directory if inter_dir isn't specified
    inter_dir = pathlib.Path(os.path.abspath(".")) if inter_dir is None else pathlib.Path(os.path.abspath(inter_dir))

    anndata_dirs, anndata_sparse = aggregate(in_dir, aggregated_filename, inter_dir)
    anndata_annot = annotate(anndata_sparse, anndata_dirs, inter_dir)
    qc_keyword_arguments = {
        "pct_count": (pct_count_lower, pct_count_upper),
        "genes_by_count": (genes_by_count_lower, genes_by_count_upper),
        "total_counts": (total_counts_lower, total_counts_upper)
    }
    anndata_filtered = quality_control(anndata_annot, inter_dir, **qc_keyword_arguments)
    postprocess(anndata_filtered, inter_dir)


if __name__ == "__main__":
    preprocess()
