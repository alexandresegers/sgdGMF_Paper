#' 
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
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
dec.sce <- readRDS(file = "Output_files/dec.sce.RDS")
rowData(tenx)$bio.var <- dec.sce$bio
hvg <- (rowData(tenx)$Ensembl[order(rowData(tenx)$bio.var, decreasing = TRUE)])[1:500]

tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg,]

#' 
#' 
## -------------------------------------------------------------------------------------------------
Y <- as.matrix(t(counts(tenx_filtered)))
X <- model.matrix(~1, colData(tenx_filtered))

#' 
## -------------------------------------------------------------------------------------------------
set.seed(100)
family = poisson()
method = "sgd"
sampling = "block"
control = list(maxiter = 10000, size = c(1000,250))

t0 <- Sys.time()
fit <- sgdGMF::sgdgmf.fit(Y, X, ncomp = 10, family = family,
                          method = method, sampling = sampling, 
                          control.alg = control,
                          control.init = list(method = "ols", type = "link"))
t1 <- Sys.time()
time <- t1 - t0  

saveRDS(object = time, file = "Output_files/CaseStudy_modelfit_full_time.RDS")
saveRDS(object = fit, file = "Output_files/CaseStudy_modelfit_full.RDS")

#' 
