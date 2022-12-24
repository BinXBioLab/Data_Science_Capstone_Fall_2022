### Note: this script is very memory intensive

# Set max virtual memory size large enough to handle in-memory operations
Sys.setenv("R_MAX_VSIZE" = 48000000000)

# Check that we have the correct number of arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop(
        "No arguments. Format: Rscript --vanilla counts metadata prefix",
        call. = FALSE
    )
} else if (length(args) == 1) {
    stop(
        "Only 1 argument. Format: Rscript --vanilla counts metadata prefix",
        call. = FALSE
    )
} else if (length(args) == 2) {
    stop(
        "Only 2 arguments. Format: Rscript --vanilla counts metadata prefix",
        call. = FALSE
    )
} else if (length(args) > 3) {
    stop(
        "Too many arguments. Format: Rscript --vanilla counts metadata prefix",
        call. = FALSE
    )
}

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
output_prefix <- as.character(args[3])

print(paste("Reading counts from:", counts_path, sep = " "))
counts <- arrow::read_parquet(counts_path)

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

print("Running Optimize ALS")
sdata <- RunOptimizeALS(sdata, k = 20, lambda = 5, split.by = "batch")

print("Quantile Normalization")
sdata <- RunQuantileNorm(sdata, split.by = "batch")

print("Writing Outputs to Disk")
arrow::write_parquet(
    x = as.data.frame(sdata@reductions$iNMF@cell.embeddings),
    sink = paste(output_prefix, ".liger.parquet", sep = "")
)

arrow::write_parquet(
    as.data.frame(sdata@reductions$iNMF_raw@cell.embeddings),
    sink = paste(output_prefix, ".liger.usage.parquet", sep = "")
)

arrow::write_parquet(
    as.data.frame(sdata@reductions$iNMF_raw@feature.loadings),
    sink = paste(output_prefix, ".liger.weights.parquet", sep = "")
)