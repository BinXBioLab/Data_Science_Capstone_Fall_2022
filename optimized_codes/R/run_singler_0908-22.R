##############################################################
# run singleR for ref cell annotation transfer to target cells
# Bin Xu 09-08-22
# To avoid a memory issue. We did not include the permutation test
# 
#############################################################
# Inputs:
# NOTE: SingleR() expects log-counts, this script expects raw counts
# Log transformations automatically done in the script

# ref = nowakowski cpm matrix in feather format generated from python
# mat = out scRNAseq data cpm matrix in feather format generated from python
# labels = cell type annotation of all nowakowski cells (The
# definition of cell types can be modified in this file) from python
# filter = 'noglyc' to remove the Glyc cells, which might be related to stress

library(SingleR)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(BiocParallel)
library (lubridate)
library(arrow)

#library(feather) # this will cause segfault problem. use arrow instead
    ## need to conda install -c conda-forge r-arrow
    ## version 8.0.0 works

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
    
refset <- "nowakowski"
label_def <- "med"
targetset <- "setd1a_b2"
rundate <- today()

## define inputs and outputs here
input1 <- "nowakowski_090422_v1.cpm.fth"
input2 <- "setd1a_b2_090422_v1.cpm.fth"
input3 <- "nowakowski.label_med.csv"

out1 <- paste(targetset, refset, label_def, rundate, "singler.csv", sep="_")
out2 <- paste(targetset, refset, label_def, rundate, "noglyc_singler.csv", sep="_")

# two runs. One with and one without the filte
for (test in c("", "noglyc")) {
    rm(list = setdiff(ls(), "test"))


# if errors, check the format of input files.
# It should use cellname as rownames and gene names as colnames
    ref <- as.data.frame(read_feather(input1))
    rownames(ref) <- ref[,"Cell"]
    ref <- ref[,-1]
    mat <- as.data.frame(read_feather(input2))
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]

    labels <- read.table(
        sep=',',
        header=TRUE,
        check.names=FALSE,
        input3,
        row.names=1,
    )
    if(test == 'noglyc') {
        labels <- subset(labels, Ref != 'Glyc')
    }

    # filter out glyc cells annotation if needed
    if (refset == 'nowakowski') {
        cells <- intersect(rownames(ref), rownames(labels))
        ref <- ref[cells,]
    }

    # Transpose and log1p of the ref matrix.
    if (refset == 'nowakowski') {
        ref <- log1p(t(ref))
    }

   # Transpose and log1p of the ref matrix.
    mat <- log1p(t(mat))
   
   # filter the labels
    if (refset == 'nowakowski') {
        labels.full <- labels[colnames(ref),]
    } 

    # run singleR to map the ref cell annotation to target cells
    pred <- SingleR(
        test=mat,
        ref=ref,
        labels=labels.full,
        de.method='wilcox',
        BPPARAM=MulticoreParam(workers=4+1),
    )

    if (test != 'noglyc') {
        write.table(
            pred,
            file=out1,
            sep=',',
            quote=FALSE,
        )
    } else {
        write.table(
            pred,
            file=out2,
            sep=',',
            quote=FALSE,
        ) 
    }

}
