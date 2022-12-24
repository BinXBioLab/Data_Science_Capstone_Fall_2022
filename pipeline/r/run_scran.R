library(scran)
library(tidyverse)
library(Matrix)

print("Parsing command line arguments")
args <- commandArgs(trailingOnly = TRUE)
format <- "Format: Rscript --vanilla input_file.mtx output_file.mtx"

if (length(args) == 0) {
    stop(
        paste("No arguments.", format, sep = " "),
        call. = FALSE
    )
} else if (length(args) == 1) {
    stop(
        paste("Too few arguments.", format, sep = " "),
        call. = FALSE
    )
} else if (length(args) > 2) {
    stop(
        paste("Too many.", format, sep = " "),
        call. = FALSE
    )
}

# Store input/output file paths
input_path <- args[1]
output_path <- args[2]

print("Creating matrix and clustering")
mat <- Matrix::readMM(input_path)
cl <- quickCluster(mat)

print("Computing sum factors")
mysce <- computeSumFactors(
    SingleCellExperiment(list(counts = mat)),
    clusters = cl
)

print("Write array of double values to disk")
sf <- sizeFactors(mysce)
write.table(sf, file = output_path, col.names = FALSE, row.names = FALSE)
