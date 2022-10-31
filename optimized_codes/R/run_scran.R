library(scran)
library(tidyverse)
library(Matrix)

mat <- Matrix::readMM("optimized_codes/scripts/test.mtx")
cl <- quickCluster(mat)

mysce <- computeSumFactors(
    SingleCellExperiment(list(counts=mat)),
    clusters = cl
)

sf <- sizeFactors(mysce)