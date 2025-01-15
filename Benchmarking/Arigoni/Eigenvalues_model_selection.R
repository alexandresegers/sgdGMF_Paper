#' 
## --------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)

#' 
## --------------------------------------------------------------------------------------------------
sce <- readRDS(file = "Data/BE1/final/sce.RDS")
sce <- sce[(rownames(sce)[order(rowData(sce)$biological.var, 
                                decreasing = TRUE)])[1:500],]
Y <- as.matrix(t(counts(sce)))


#' 
#' 
## --------------------------------------------------------------------------------------------------
set.seed(100)
eigs <- sgdgmf.rank(Y,
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

saveRDS(eigs, file = "Output_models/eigenvalues.RDS")

#' 
#' 
#' 
