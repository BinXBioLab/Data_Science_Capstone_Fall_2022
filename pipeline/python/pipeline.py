from __future__ import annotations
import pdb
import sys
import pickle
import scanpy as sc
import git
import pathlib
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from helper import write_anndata, create_color_palette
import wget
from scipy.sparse import csr_matrix
from tqdm import tqdm
import gseapy
import click
from joblib import parallel_backend

sns.set_style('white')
plt.rcParams['savefig.facecolor'] = 'w'
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.set_figure_params(
    scanpy=True, 
    dpi=120, 
    dpi_save=150, 
    frameon=True, 
    vector_friendly=True, 
    fontsize=14
)
sc.settings.figdir = os.path.join(outDir, 'fig_supp')
sc.settings.n_jobs = 4

@click.group()
def cli():
    """
    scRNA Analysis Pipeline: Automated sequence of data analyses for scRNA sequencing data
    
    Visit https://github.com/BinXBioLab/Data_Science_Capstone_Fall_2022 for usage instructions and help
    """

def generate_figures(anndata: sc.AnnData) -> None:
    # Compute neighborhood graph
    # use_rep means use X_liger result to represent X for calculation
    # other adjustable parameters that control uMAP output are n_neighbors, metrics(distance), min_dist
    with parallel_backend('threading', n_jobs=4):
        sc.pp.neighbors(
            anndata,
            n_neighbors=15,
            random_state=42,
            use_rep= 'X_liger'
        )

    # conduct umap
    with parallel_backend('threading', n_jobs=4):
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

def rank_genes_and_save(adata: sc.AnnData) -> None:
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', show=False)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
    write_anndata(adata, interDir, "pasca_log1p_liger.h5ad")

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
    meta = pd.read_csv(meta_filepath, sep = '\t')
    h_ad.obs = meta
    anndata = sc.AnnData(h_ad)

    return anndata

def visualize_scrna(adata):
    markers = ['DGCR8', 'RANBP1', 'DGCR2', 'COMT', 'CRKL']
    sc.pl.dotplot(
        adata, 
        markers, 
        groupby='sampleID', 
        dendrogram=False, 
        save='del_gene.pdf',
        show=False
    )

    # Plot UMAP based on certain genes
    sc.pl.umap(
        adata,
        color=['leiden','PAX6','SOX2','FOXG1','NKX2-1'],
        cmap='viridis',
        size=5,
        wspace=.3,
        save='key_genes.pdf',
        show=False
    )

    # Plot UMAP based on nowakowski_noglyc_coarse
    sc.pl.umap(
        adata,
        color=['nowakowski_noglyc_coarse' ],
        cmap='viridis',
        size=5,
        wspace=.3,
        save='.Annotation_coarse.pdf',
        show=False
    )

    # Not sure where my_colors is used...
    with open('/content/drive/My Drive/Colab Notebooks/my_palette.txt', 'rb') as f:
        my_colors = pickle.load(f)

    for i in tqdm(adata.obs.condition.cat.categories):
        sc.set_figure_params(scanpy=True, dpi=200, dpi_save=300, frameon=True, vector_friendly=True, fontsize=14)
        sc.pl.umap(adata[adata.obs.condition==i,:],color=['leiden'], legend_fontsize=6,legend_loc="on data",legend_fontoutline=1)
        sc.pl.umap(adata[adata.obs.condition==i,:],color=['leiden'], legend_fontsize=6,) 

    for i in tqdm(adata.obs.batch.cat.categories):
        sc.set_figure_params(scanpy=True, dpi=200, dpi_save=300, frameon=True, vector_friendly=True, fontsize=14) 
        sc.pl.umap(adata[adata.obs.batch==i,:],color=['NEUROD2'],legend_fontsize=5, palette=my_colors, legend_loc="on data",legend_fontoutline=1, title=i)

    for i in tqdm(adata.obs.condition.cat.categories):
        sc.set_figure_params(scanpy=True, dpi=250, dpi_save=600, frameon=True, vector_friendly=True, fontsize=14) 
        sc.pl.umap(adata[adata.obs.condition==i,:],color=['nowakowski_noglyc_coarse'],legend_fontsize=7, palette=my_colors, legend_loc="on data",legend_fontoutline=1, title=i, save=f'_{i}_nowakowski_noglyc_coarse.png')

    for i in tqdm(adata.obs.condition.cat.categories):
        sc.set_figure_params(scanpy=True, dpi=200, dpi_save=150, frameon=True, vector_friendly=True, fontsize=14) 
        sc.pl.umap(adata[adata.obs.condition==i,:],color=['SETD1A'],legend_fontsize=5, palette=my_colors, legend_loc="on data",legend_fontoutline=1, title=i, save=f'.celltypes.pdf',)
    
def quantitative_cell_composition(adata, gene="SETD1A", palette=None):
    sc.set_figure_params(
        scanpy=True, 
        dpi=300, 
        dpi_save=600, 
        frameon=True, 
        vector_friendly=True, 
        fontsize=7
    )

    for i in adata.obs.condition.cat.categories:
        sc.pl.violin(
            adata[adata.obs.condition==i,:],
            keys= gene, 
            groupby="nowakowski_noglyc_coarse", 
            palette=palette,
            rotation=90,
            size=0.5,
            ylabel=f'{gene} Expression level',
            save=f'_{gene}_{i}_celltype.png'
        )

    for i in adata.obs.condition.cat.categories:
        sc.set_figure_params(scanpy=True, dpi=120, dpi_save=600, frameon=True, vector_friendly=True, fontsize=14) 
        sc.pl.umap(
            adata[adata.obs.condition==i,:],
            color=['nowakowski'], 
            palette=palette, 
            title=i, 
            save=f'{i}_nowakowski.pdf',
            legend_loc="on data", 
            legend_fontsize=8,
            legend_fontoutline=1,
        )

    for i in adata.obs.condition.cat.categories:
        sc.pl.dotplot(
            adata[adata.obs.condition== i,:],
            var_names=gene,
            groupby='nowakowski',
            title=i, 
            save=f'{i}_nowakowski.pdf'
        )

    for i in adata.obs.condition.cat.categories:
        sc.pl.dotplot(
            adata[adata.obs.condition== i,:],
            var_names= gene,
            groupby='leiden',
            title=i, 
            save=f'{i}_leiden.pdf'
        )

# TODO: Turn print statements into pandas DataFrame and save it
def rank_gene_analysis(adata: sc.AnnData, gene='SETD1A', count_min=1) -> sc.AnnData:
    # run the rank gene analysis
    sc.tl.rank_genes_groups(adata, 'leiden')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    new_data = {
        group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names','logfoldchanges','pvals_adj']
    }
    DE_genes_tmp = pd.DataFrame(new_data)
    DE_genes_tmp.to_parquet(os.path.join(interDir, f'{prefix}.log1p_liger_singleR_Rank_leiden.parquet'))

    # counts per genotype
    for i in adata.obs.condition.cat.categories:
        for j in adata.obs.sampleID.cat.categories:
            if not (adata[(adata.obs.condition==i) & (adata.obs.sampleID==j),:].obs['nowakowski_med'].value_counts().empty):
                c = adata[(adata.obs.condition==i) & (adata.obs.sampleID==j),:].obs['nowakowski_med'].value_counts()
                print ("<<" + j)
                print (c)

    # Get cells with high targetGene expression X > mycut
    # generate a filter in the obs group
    adata.obs['targetGene'] = (adata[:,[gene]].X > count_min).sum(1) 
    
    # subset selection
    adata_test = adata[adata.obs['targetGene'] > 0]
    
    # recompute number of genes expressed per cell
    adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1)

    # counts per genotype
    for i in adata_test.obs.condition.cat.categories:
        print (">>" + gene + "_" + i)
        c = adata_test[adata_test.obs.condition==i,:].obs['nowakowski_noglyc_coarse'].value_counts()
        print (c)

    return adata

def mast_preprocessing(adata, cutoff=0, gene = 'SETD1A'):
    # backup new X data
    adata.layers['bck_X'] = adata.X

    # read back the raw X data
    adata.X = adata.layers['counts']

    # normalized and log2
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=2)

    # subset the adata for the cell type we are interested in
    cellType = ['EN']
    adata_test = adata[adata.obs['nowakowski'].isin(cellType)]
    adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1) # recompute number of genes expressed per cell

    # get the cells with high targetGene expression X > cutoff
    # generate a filter in the obs group
    adata_test.obs['targetGene'] = (adata_test[:,[gene]].X > cutoff).sum(1) 

    adata_test2 = adata_test[adata_test.obs['targetGene'] >= 1]
    
    # Ranking by t-test
    sc.tl.rank_genes_groups(adata_test2, 'condition', method='t-test', key_added = "t-test")
    sc.pl.rank_genes_groups(adata_test2, n_genes=25, sharey=False, key = "t-test")

    # Ranking by Wilcoxon
    sc.tl.rank_genes_groups(adata_test2, 'condition', method='wilcoxon', key_added = "wilcoxon")
    sc.pl.rank_genes_groups(adata_test2, n_genes=25, sharey=False, key="wilcoxon")

    # Ranking by logistic regression
    sc.tl.rank_genes_groups(adata_test2, 'condition', method='logreg',key_added = "logreg")
    sc.pl.rank_genes_groups(adata_test2, n_genes=25, sharey=False, key = "logreg")
    
    sc.pl.rank_genes_groups_dotplot(adata_test2, n_genes=10, key="logreg", groupby='condition',groups=['SP'])
    sc.pl.rank_genes_groups_stacked_violin(adata_test2, n_genes=5, key="logreg", groupby="condition")
    sc.pl.rank_genes_groups_violin(adata_test2, n_genes=10, key="logreg")

    gene_set_names = gseapy.get_library_name(database='Human')
    res = sc.get.rank_genes_groups_df(adata_test2, group=['CT','FS','SP'], key='wilcoxon',pval_cutoff=0.1)
    res.to_parquet(os.path.join(interDir, f'{prefix}.DEG_rank_genes.parquet'))

def differential_expression_analysis(adata, cutoff=0, gene = 'SETD1A'):
    return run_mast(adata)


def find_highly_variable_genes(adata):
    ad = adata # remove this later
    sc.pp.highly_variable_genes(ad, min_mean= 0.2,max_mean=10, min_disp=0.1)
    sc.pl.highly_variable_genes(ad)

    # form a new object with highly variable genes only. Assume these genes are the key genes that can explain various differences in the data
    ad_hv = ad[:, ad.var.highly_variable]
    filename = "cntnap2-s26_norm_log_hv.h5ad"
    output_path = os.path.join(interDir, filename)
    ad_hv.write(output_path)

def visualize_highly_variable_genes(adata):
    ad_hv = adata # remove this later

    # adjust or readjust the pca and umap
    sc.tl.pca(ad_hv, svd_solver='arpack')
    sc.pp.neighbors(
        ad_hv,random_state=42,
        use_rep= 'X_liger_k25'
    )
    sc.tl.umap(ad_hv,spread=3)

    # TODO: Make this more general/modular 
    sc.pl.umap(
        ad_hv, 
        color=['CNTNAP2','PAX6','NEUROD6','NEUROD2'], 
        palette='hsv', 
        legend_fontsize=8, 
        legend_loc="on data"
    )

    # Format later
    sc.pl.umap(ad_hv[ad_hv.obs['condition'] =='WT',:],color=['CNTNAP2','PAX6','NEUROD6','NEUROD2'], palette='hsv', legend_fontsize=6,title='WT : Leiden', legend_loc="on data")
    sc.pl.umap(ad_hv[ad_hv.obs['condition'] =='MUT',:],color=['CNTNAP2','PAX6','NEUROD6','NEUROD2'], palette='hsv', legend_fontsize=6,title='MUT : Leiden',legend_loc="on data")

    filename = "cntnap2-s26_norm_log_hv_leiden.h5ad"
    output_path = os.path.join(interDir, filename)
    ad_hv.write(output_path)

def human_pfc_data(
    outfile: str | pathlib.Path, 
    exprMatrix:str|None = None, 
    meta:str|None = None) -> sc.AnnData:
    
    # Download reference material
    if exprMatrix == None:
        matrix_filepath = os.path.join(interDir, "exprMatrix.tsv.gz")
        wget.download("https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz", matrix_filepath)

    if meta == None:
        meta_filepath = os.path.join(interDir, "meta.tsv")
        wget.download("https://cells.ucsc.edu/cortex-dev/meta.tsv", meta_filepath)
    
    meta = pd.read_csv(meta_filepath, sep='\t')
    exprMatrix = sc.read_text(exprMatrix, sep='\t')

    h_ad = exprMatrix.transpose()
    h_ad.obs = meta
    output_path = os.path.join(interDir, outfile)
    h_ad.write(output_path)

    return h_ad

@cli.command()
@click.option('--rscript', type=click.Path(exists=True), required=True, default=None, help='R script that will run Liger')
@click.option('--counts', type=click.Path(exists=True), required=True, default=None, help='Path to counts from preprocessing step')
@click.option('--metadata', type=click.Path(exists=True), required=True, default=None, help='Path to metadata from preprocessing step')
@click.option('--output', type=str, required=True, default=None, help='Filename for output of Liger analysis')
def run_liger(rscript, counts, metadata, output):
    click.echo("Preprocessing --> *You are here* --> post_liger_pre_singler --> ...")
    
    # Construct command line arguments to run Liger script and run it
    args = f"{rscript} {counts} {metadata} {output}"
    os.system(f"Rscript --vanilla {args}")
    
    return output

@cli.command()
@click.option('--anndata', type=click.Path(exists=True), default=None, help='h5ad file from end of preprocessing')
@click.option('--liger', type=click.Path(exists=True), default=None, help='Liger output')
@click.option('--output', type=click.Path(exists=True), default=None, help='Output file path')
def post_liger_pre_singler(anndata, liger, output):
    click.echo("Preprocessing --> Liger --> *You are here* --> Singler")
    inter = os.path.abspath(input) if input is not None else interDir
    print(f"Using input directory: {inter}")
    
    click.echo("Reading output of Liger + Preprocessing")
    df_liger = pd.read_csv(liger, index_col=0)
    adata = sc.read(anndata)

    # Generating figures based on Liger output
    click.echo("Generating figures based on Liger output")
    adata.obsm['X_liger'] = df_liger.loc[adata.obs_names, :].values
    adata = generate_figures(adata)
    
    adata.uns['log1p']["base"] = None # To avoid KeyError: 'base' with rank_gene_groups()
    rank_genes_and_save(adata)

    click.echo("Reading and processing post quality control data")
    postqc_path = os.path.join(interDir, "pasca_postqc.h5ad")
    adata = sc.read(postqc_path)
    df = pd.DataFrame.sparse.from_spmatrix(adata.X, index=adata.obs_names, columns=adata.var_names)

    # filter genes?
    click.echo("Filtering genes")
    df = df.loc[:, ((df > 0).sum(axis=0) >= 1)]

    click.echo("Downloading and processing reference data")
    ad_nowakowski = download_and_process_reference()
    sc.pp.filter_cells(ad_nowakowski, min_genes=200)
    sc.pp.filter_genes(ad_nowakowski, min_cells=1)
    genes = df.columns.intersection(ad_nowakowski.var_names)

    df = df.loc[:, genes]
    ad_nowakowski = ad_nowakowski[:, genes].copy()

    # df is the target dataframe to be exported to R 
    # generate the ref dataframe with CellID as index to be exported to R
    click.echo("Generating reference dataframe to be exported to R")
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

    click.echo("Start of part 7d in original notebook")
    ## 7d START ## 
    # Presumbly df refers to earlier df?
    res1 = (df * 1e6).divide(df.sum(axis=1), axis='rows')
    c_res1 = csr_matrix(res1)
    df1 = pd.DataFrame(c_res1.toarray(), index=df.index, columns=df.columns)
    df1.reset_index().to_parquet(output)

    res2 = (df_nowakowski * 1e6).divide(df_nowakowski.sum(axis=1), axis='rows')
    c_res2 = csr_matrix(res2)
    df2 = pd.DataFrame(c_res2.toarray(),index=df_nowakowski.index,columns=df_nowakowski.columns)
    df2.reset_index().to_parquet(os.path.join(interDir, "nowakowski_pasca.cpm.parquet")) 
    
    # Return path to input file for SingleR
    return output

@cli.command()
@click.option('--rscript', type=click.Path(exists=True), required=True, default=None, help='R script that will run SingleR')
@click.option('--output', type=click.Path(exists=True), default=None, required=True, help='Output file path')
def run_singler(rscript, input_file, output):
    click.echo("post_liger_pre_singler --> *You are here* --> post_singler_pre_mast --> ...")
    args = f"{rscript} {input_file} {output}"
    cmd = f"Rscript --vanilla {args}"
    os.command(cmd)
    
    return output

@cli.command()
@click.option('--anndata', type=click.Path(exists=True), default=None, help='h5ad file from previous step')
@click.option('--annotation', type=click.Path(exists=True), default=None, help='Annotation file to augment anndata')
@click.option('--inter', type=click.Path(exists=True), default=None, help='Output file path')
def post_singler_pre_mast(anndata, annotation, inter):
    click.echo("... --> Singler --> *You are here* --> MAST")
    
    ## 7f START ##
    adata = sc.read(anndata)
    annotation = pd.read_csv(annotation, index_col=0)
    adata.obs[f'nowakowski_med'] = annotation.loc[adata.obs_names, 'labels'] 
    del adata.obs['nowakowski_noglc'] # not sure why we're deleting this
    ## 7f END ##

    ## Part 6 START ##
    adata.write(os.path.join(inter, "pasca.log1p_liger_med_singleR_noglyc.h5ad"))
    adata.write(os.path.join(inter, "tangram_gpu_input.h5ad"))
    adata.write(os.path.join(inter, "cellphonedb_input.h5ad"))
    ## Part 6 END ##

    ## Part 7 START ##
    visualize_scrna(adata)
    ## Part 7 END ##

    ## Part 8 START ##
    quantitative_cell_composition(adata)
    ## Part 8 END ##

    ## Part 9 START ##
    adata = rank_gene_analysis(adata)
    ## Part 9 END ##
    
    adata = mast_preprocessing(adata)
    
    return adata

@cli.command()
@click.option('--rscript', type=click.Path(exists=True), required=True, help='R script that will run MAST')
@click.option('--anndata', type=click.Path(exists=True), required=True, help='Path to anndata')
@click.option('--metadata', type=click.Path(), required=True, help='Filepath to store intermediate metadata for MAST')
@click.option('--output', type=click.Path(), default=None, required=True, help='Output file path')
def run_mast(rscript, anndata, metadata, output) -> None:
    # TODO: Turn cell_types into a CLI argument and function argument
    cell_types=['EN-PFC1','EN-PFC2','EN-PFC3','EN-V1-1','EN-V1-2','EN-V1-3']
    adata = sc.read(anndata)
   
    # subset the adata for the cell type we are interested in
    adata_test = adata[adata.obs['nowakowski'].isin(cell_types)]
    adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1) # recompute number of genes expressed per cell

    # produce mast data csv
    adata_test.to_df().to_csv(metadata)

    # produce mast metadata csv
    #, ... additional covariates present in obs slot of adata 
    # that you may want to use ...  ]
    keys = ['n_genes', 'condition', 'batch'] 
    adata_test.obs[keys].to_csv(metadata)
    
    cmd = f"Rscript --vanilla {rscript} {metadata} {output}"
    os.system(cmd)

    ent_de = pd.read_csv(output)
    return ent_de

@cli.command()
@click.option('--ent_de', type=str, default=None, help='Output from differential gene analysis via MAST R package')
def post_mast(ent_de):
    click.echo("--> post_singler_pre_mast --> MAST --> *You are here*")
    ent_de = pd.read_csv(ent_de)
    ent_de.index = ent_de.primerid

    ## import annotation from the biomart database
    annot = sc.queries.biomart_annotations(
                "hsapiens",
                ['hgnc_symbol',"chromosome_name","start_position", "end_position"],
            ).set_index("hgnc_symbol")

    ent_de_ann = ent_de.merge(annot, left_index=True, right_index=True)
    ent_de_ann.sort_values(by=['FDR'], inplace=True)

    filename = f"{prefix}.log1p_liger_singleR_nowakowski_3cellTypes_MAST.parquet"
    output_path = os.path.join(interDir, filename)
    ent_de_ann.to_parquet(output_path)
    return

if __name__ == "__main__":
    git_repo = git.Repo(".", search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")

    topDir = os.path.join(git_root, "raw_data/Pasca_scRNAseq/")
    prefix = 'pasca'
    outDir = pathlib.Path(topDir,"outputs")
    interDir = pathlib.Path(topDir, "inter-test")
    
    cli()
