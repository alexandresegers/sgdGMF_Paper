#' 
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)

#' 
## -------------------------------------------------------------------------------------------------
tenx <- loadHDF5SummarizedExperiment(
                             dir =  "TENxBrainData_preprocessed_default", 
                             prefix="")

#' 
#' 
## -------------------------------------------------------------------------------------------------
dec.sce <- readRDS(file = "Output_files/dec.sce.RDS")
rowData(tenx)$bio.var <- dec.sce$bio
hvg <- (rowData(tenx)$Ensembl[order(rowData(tenx)$bio.var, decreasing = TRUE)])[1:500]

tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg,]

#' 
#' 
#' 
## -------------------------------------------------------------------------------------------------
Y <- as.matrix(t(counts(tenx_filtered)))

#' 
#' 
## -------------------------------------------------------------------------------------------------
ranks <- sgdgmf.rank(Y,
    X = NULL,
    Z = NULL,
    maxcomp = 50,
    family = poisson(),
    method = c("onatski"),
    type.reg = c("ols"),
    type.res = c("link"),
    maxiter = 10,
    parallel = FALSE,
    nthreads = 1)

saveRDS(ranks, file = "Output_files/eigenvalues_sgdGMF_link.RDS")

