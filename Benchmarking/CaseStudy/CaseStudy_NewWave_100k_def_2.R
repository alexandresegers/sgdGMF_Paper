#' 
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)
library(dplyr)
library(NewWave)

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
samples <- sample(x = colnames(tenx), size = 100000, replace = FALSE)
tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg, samples]

#' 
#' 
## -------------------------------------------------------------------------------------------------
Y <- as.matrix(t(counts(tenx_filtered)))
X <- model.matrix(~1, colData(tenx_filtered))

#' 
## -------------------------------------------------------------------------------------------------
## NBWaVE
fit.nbwave = function (
    y, x = NULL, z = NULL, ncomp = 10, family = poisson(),
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  idx = which(apply(y, 2, sd) == 0)
  if (length(idx) > 0) {
    for (j in idx) {
      y[sample(1:n, 1),j] = 1
    }
  }

  # NBWaVE model fitting
  time0 = proc.time()
  fit = NewWave::newFit(
    Y = t(y),
    X = x,
    V = z,
    K = ncomp,
    commondispersion = TRUE,
    children = 4,
    verbose = verbose) %>%
    suppressWarnings()
  timef = proc.time()



  # Output
  list(
    model = "NBWaVE",
    fit = fit,
    time = timef - time0)
}

#' 
## -------------------------------------------------------------------------------------------------
set.seed(1002)

ncomp <- 10

fit = fit.nbwave(y = Y, x = X, ncomp = ncomp)

saveRDS(object = fit, file = "Output_files/NewWave_100k_def_2.RDS")

#' 
