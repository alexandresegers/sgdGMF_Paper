#' 
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)
library(dplyr)
library(sgdGMF)

#' 
#' 
## -------------------------------------------------------------------------------------------------
tenx <- loadHDF5SummarizedExperiment(
                             dir =  "TENxBrainData_preprocessed_default", 
                             prefix="")

#' 
#' 
## -------------------------------------------------------------------------------------------------
set.seed(100)
dec.sce <- readRDS(file = "Output_files/dec.sce.RDS")
rowData(tenx)$bio.var <- dec.sce$bio
hvg <- (rowData(tenx)$Ensembl[order(rowData(tenx)$bio.var, decreasing = TRUE)])[1:500]
samples <- sample(x = colnames(tenx), size = 300000, replace = FALSE)
tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg, samples]

#' 
#' 
## -------------------------------------------------------------------------------------------------
Y <- as.matrix(t(counts(tenx_filtered)))
X <- model.matrix(~1, colData(tenx_filtered))

#' 
#' 
## -------------------------------------------------------------------------------------------------
family = poisson()
ncomp <- 10


for(i in c(1:5)){
  t0 <- Sys.time()
  fit <- sgdGMF::sgdgmf.fit(Y=Y, X=X, ncomp=ncomp, family=family, method = "sgd",
                            sampling = "block",
                            control.alg = list(maxiter = 3000, size = c(1000,250)),
                            control.init = list(method = "ols", type = "link"))
  t1 <- Sys.time()
  time <- t1 - t0  
  saveRDS(object = list("fit" = fit, "time" = time), 
        file = paste0("Output_files/sgdgmffit_300k_samples_3000_maxiter_1000_250_size_",i ,".RDS"))
  rm(fit)

}


#' 
