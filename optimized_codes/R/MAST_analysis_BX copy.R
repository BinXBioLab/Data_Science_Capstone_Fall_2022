# Check that we have the correct number of arguments
if (length(args)==0) {
  stop("No arguments found, you must pass 1 file for the script", call.=FALSE)
} else if (length(args)>1) {
  stop("You must pass only 1 file for the script", call.=FALSE)
}

library(mast)
library(arrow)
library(anndata)
library(SingleCellExperiment)

# Get path of anndata file from command line
anndata_path = args[1]
adata = anndata::read_h5ad(anndata_path)

# Get the data from the anndata object into a SingleCellExperiment object
sca <- SceToSingleCellAssay(adata, class = 'SingleCellAssay')

# TODO: Find some way to pass in conditions in command line arguments
conditions <- c("CT", "FS", "targetGene1")

# Loop over all the conditions and perform differential gene analysis on each condition
for (condition in conditions) {
  # Set the reference set
    cond <- factor(colData(sca)$condition)
    cond<-relevel(cond, condition)
    colData(sca)$condition<-cond

    sca_ent_filt = sca[rowSums(assay(sca)) != 0, ]

    ## This the most time consuming steps
    zlmCond_ent <- zlm(formula = ~condition + n_genes, sca=sca_ent_filt)

    ## This the second most time consuming steps
    summaryCond_ent <- summary(zlmCond_ent, doLRT=paste('condition', condition, sep=''))
    summaryDt_ent <- summaryCond_ent$datatable

    ## get results
    result_ent <- merge(
        summaryDt_ent[contrast=='conditionFS' & component=='H',
        .(primerid, `Pr(>Chisq)`)], #P-vals
        summaryDt_ent[contrast=='conditionFS' & component=='logFC', 
        .(primerid, coef, ci.hi, ci.lo)],
        by='primerid'
    ) #logFC coefficients

    #Correct for multiple testing (FDR correction) and filtering
    result_ent[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    ent_de = result_ent[result_ent$FDR <= 1,, drop=F]
    ent_de = ent_de[order(ent_de$FDR),]

    # Save the ent_de table
    anndata::write_h5ad(ent_de, paste(condition, "_ent_de.h5ad", sep=""))
}