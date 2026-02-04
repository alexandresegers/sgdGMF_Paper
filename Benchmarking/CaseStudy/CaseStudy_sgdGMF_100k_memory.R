#'
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)
library(dplyr)
library(sgdGMF)
library(peakRAM)

args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[[1]])
print(i)

#'
#'
## -------------------------------------------------------------------------------------------------
tenx <- loadHDF5SummarizedExperiment(
  dir =  "../TENxBrainData_preprocessed_default",
  prefix="")

#'
#'
## -------------------------------------------------------------------------------------------------
set.seed(100)
dec.sce <- readRDS(file = "../Output_files/dec.sce.RDS")
rowData(tenx)$bio.var <- dec.sce$bio
hvg <- (rowData(tenx)$Ensembl[order(rowData(tenx)$bio.var, decreasing = TRUE)])[1:500]
samples <- sample(x = colnames(tenx), size = 100000, replace = FALSE)
tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg, samples]

#'
#'
## -------------------------------------------------------------------------------------------------
#'
#'
## -------------------------------------------------------------------------------------------------
family = poisson()
ncomp <- 10
Y <- as.matrix(t(counts(tenx_filtered)))
X <- model.matrix(~1, colData(tenx_filtered))


set.seed(100+i)

t0 <- Sys.time()
memory = peakRAM::peakRAM(
  fit <- sgdGMF::sgdgmf.fit(Y = Y,
                            X=X, ncomp=ncomp, family=family, method = "sgd",
                          sampling = "block",
                          control.alg = list(maxiter = 1000, size = c(1000,250),
                                             savedata = FALSE, normalize = FALSE),
                          control.init = list(method = "light", type = "link",
                                              normalize = FALSE))
  )
t1 <- Sys.time()
time <- t1 - t0
saveRDS(object = list("fit" = fit,
                      "time" = time,
                      "memory" = memory),
        file = paste0("Output_files/sgdgmffit_100k_memory",i ,".RDS"))






#'
