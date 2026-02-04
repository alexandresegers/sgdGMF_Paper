#'
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)
library(dplyr)
library(COAP)
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
samples <- sample(x = colnames(tenx), size = 75000, replace = FALSE)
tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg, samples]

#'
#'
## -------------------------------------------------------------------------------------------------
Y <- as.matrix(t(counts(tenx_filtered)))
X <- model.matrix(~1, colData(tenx_filtered))

#'
## -------------------------------------------------------------------------------------------------
## COAP
fit.coap = function(
    y, x = NULL, z = NULL, ncomp = 2,
    maxiter = 1000, tol = 1e-5,
    verbose = FALSE, ncores = 4) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # COAP model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- COAP::RR_COAP(
      X_count = y,
      Z = x,
      q = ncomp,
      epsELBO = tol,
      maxIter = maxiter,
      verbose = verbose,
      joint_opt_beta = FALSE,
      fast_svd = FALSE) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Output
  list(
    model = "COAP",
    fit = fit,
    time = timef - time0,
    memory = memory)
}


#'
## -------------------------------------------------------------------------------------------------
ncomp <- 10

set.seed(100+i)

fit <- fit.coap(
  y = Y, x = X, ncomp = ncomp,
  maxiter = 1000, tol = 1e-5,
  verbose = FALSE, ncores = 4)

saveRDS(object = fit, file = paste0("Output_files/COAP_tol5_75k_def_", i, ".RDS"))
