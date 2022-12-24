library(scran)
library(tidyverse)
library(Matrix)

# Check that we have the correct number of arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop(
        "No arguments. Format: Rscript --vanilla counts metadata prefix",
        call. = FALSE
    )
} else if (length(args) > 1) {
    stop(
        "Too many arguments. Format: Rscript --vanilla counts metadata prefix",
        call. = FALSE
    )
}

mat <- Matrix::readMM("optimized_codes/scripts/test.mtx")
cl <- quickCluster(mat)

mysce <- computeSumFactors(
    SingleCellExperiment(list(counts = mat)),
    clusters = cl
)

sf <- sizeFactors(mysce)
