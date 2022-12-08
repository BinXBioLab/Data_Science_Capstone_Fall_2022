### For Liger batch correction, I will run this in R interactively!
### It will consume ~40Gb of memory in peak, so it has to been run on a server with >40GB memeory
# The first step is to setup a Renv with many packages. There are 3 ways to try installatoin different packages
# check "setup_liger.R" for installation of all needed packages and env settings

# at command line run
# source ~/.bashrc

# conda activate r_env

# check if the input file names in this Rscript are correct below

# Rscript run_liger.R

# or for interactively 

# R

### Francesco script start here
# Set max virtual memory size large enough to handle in-memory operations
Sys.setenv('R_MAX_VSIZE'=32000000000)

# Check that we have the correct number of arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop(
        "No arguments found. Format: Rscript --vanilla <counts> <metadata>",
        call. = FALSE
    )
} else if (length(args) == 1) {
    stop(
        "Only 1 argument found. Format: Rscript --vanilla <counts> <metadata>",
        call. = FALSE
    )
} else if (length(args) > 2) {
    stop(
        "Too many arguments. Format: Rscript --vanilla <counts> <metadata>",
        call. = FALSE
    )
}

# rm(list = ls())
library(patchwork)
library(rliger)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)
library(arrow)
theme_set(theme_cowplot())

counts_path <- as.character(args[1])
metadata_path <- as.character(args[2])

# Two inputs from the Python pipeline
# Note this one may not work
# data <- read.table(
#     sep = ",",
#     header = TRUE,
#     check.names = FALSE,
#     counts_path,
#     row.names = 1,
# )

print(paste("Reading counts from:", counts_path, sep = " "))
counts <- arrow::read_parquet(counts_path)
# Garbage collect counts
gc()
counts_transpose <- t(counts)

print(paste("Reading metadata from:", metadata_path, sep = " "))
metadata <- read.table(
    sep = ",",
    header = TRUE,
    check.names = FALSE,
    metadata_path,
    row.names = 1,
)

print("Creating Seurat Object")
sdata <- CreateSeuratObject(
    counts_transpose,
    project = "FULL",
    assay = "RNA",
    meta.data = metadata,
)

print("Normalizing Data")
sdata <- NormalizeData(sdata)

print("Finding Variable Features")
sdata <- FindVariableFeatures(sdata)

print("Scaling Data")
sdata <- ScaleData(sdata, split.by = "batch", do.center = FALSE)

# run with default parameter k=20 and lambda=5.
# You might want to explore the parameters' values.
# In the original analysis I used the default,
# but probably you can get more information out by tuning them.
print("Running Optimize ALS")
sdata <- RunOptimizeALS(sdata, k = 20, lambda = 5, split.by = 'batch')

# deprecated
#sdata <- RunQuantileAlignSNF(sdata, split.by='batch')

print("Quantile Normalization")
sdata <- RunQuantileNorm(sdata, split.by = 'batch')

#saveRDS(sdata, file='setd1a_b4.liger.rds')
#sdata <- readRDS('setd1a_b4.liger.rds')

### only this file will be imported into adata
prefix <- "pasca"

print("Writing Outputs to Disk")
arrow::write_parquet(
    x = as.data.frame(sdata@reductions$iNMF@cell.embeddings),
    sink = paste(prefix, ".liger.parquet", sep = "")
)

arrow::write_parquet(
    as.data.frame(sdata@reductions$iNMF_raw@cell.embeddings),
    sink = paste(prefix, ".liger.usage.parquet", sep = "")
)

arrow::write_parquet(
    as.data.frame(sdata@reductions$iNMF_raw@feature.loadings),
    sink = paste(prefix, ".liger.weights.parquet", sep = "")
)

# write.table(
#     as.data.frame(sdata@reductions$iNMF@cell.embeddings),
#     file=paste(prefix, ".liger.csv", sep=''),
#     sep=",",
#     quote=FALSE,
# )

# ## extra files saved.
# write.table(
#     as.data.frame(sdata@reductions$iNMF_raw@cell.embeddings),
#     file=paste(prefix, ".liger.usage.csv", sep=""),
#     sep=",",
#     quote=FALSE,
# )

# write.table(
#     as.data.frame(sdata@reductions$iNMF_raw@feature.loadings),
#     file=paste(prefix, ".liger.weights.csv", sep=""),
#     sep=",",
#     quote=FALSE,
# )