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
rm(list=ls())

library(rliger)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Two inputs from the Python pipeline
data <- read.table(
    sep=',',
    header=TRUE,
    check.names=FALSE,
    'pasca_log1p.csv',
    row.names=1,
)

meta.data <- read.table(
    sep=',',
    header=TRUE,
    check.names=FALSE,
    'pasca_log1p.metadata.csv',
    row.names=1,
)

sdata <- CreateSeuratObject(
    t(data),
    project='FULL',
    assay='RNA',
    meta.data=meta.data,
)

sdata <- NormalizeData(sdata)
sdata <- FindVariableFeatures(sdata)
sdata <- ScaleData(sdata, split.by='batch', do.center=FALSE)
# run with default parameter k=20 and lambda=5. You might want to explore the parameters' values. In the original analysis I used the default, but probably you can get more information out by tuning them.
sdata <- RunOptimizeALS(sdata, k=20, lambda=5, split.by='batch')
# deprecated
#sdata <- RunQuantileAlignSNF(sdata, split.by='batch')
sdata <- RunQuantileNorm(sdata,split.by='batch')

#saveRDS(sdata, file='setd1a_b4.liger.rds')
#sdata <- readRDS('setd1a_b4.liger.rds')


### only this file will be imported into adata
prefix <- 'pasca'

write.table(
    as.data.frame(sdata@reductions$iNMF@cell.embeddings),
    file=paste(prefix, '.liger.csv', sep=''),
    sep=',',
    quote=FALSE,
)

## extra files saved.
write.table(
    as.data.frame(sdata@reductions$iNMF_raw@cell.embeddings),
    file=paste(prefix, '.liger.usage.csv', sep=''),
    sep=',',
    quote=FALSE,
)

write.table(
    as.data.frame(sdata@reductions$iNMF_raw@feature.loadings),
    file=paste(prefix, '.liger.weights.csv', sep=''),
    sep=',',
    quote=FALSE,
)

