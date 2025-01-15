#' 
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)
library(sgdGMF)
library(dplyr)
library(glmpca)

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
samples <- sample(x = colnames(tenx), size = 200000, replace = FALSE)
tenx_filtered <- tenx[rowData(tenx)$Ensembl %in% hvg, samples]

#' 
#' 
## -------------------------------------------------------------------------------------------------
Y <- as.matrix((counts(tenx_filtered)))
X <- model.matrix(~1, colData(tenx_filtered))

#' 
## -------------------------------------------------------------------------------------------------
fit.glmpca = function (
    y, x = NULL, z = NULL, ncomp = 10, family = poisson(),
    method = c("avagrad", "fisher"),
    verbose = FALSE, train = NULL, test = NULL) {

  m = nrow(y)
  n = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # glmPCA model fitting
  time0 = proc.time()
  fit = glmpca::glmpca(
    Y = y, X = x, Z = z,
    L = ncomp,
    fam = "poi",
    optimizer = method,
    ctl = list(verbose = verbose)) %>%
    suppressWarnings()
  timef = proc.time()


 

  # Output
  list(
    model = switch(method, "avagrad" = "AvaGrad", "fisher" = "Fisher"),
    fit = fit,
    time = timef - time0)
}


#' 
## -------------------------------------------------------------------------------------------------
family = poisson()
ncomp <- 10

for(i in c(1:5)){
  fit <- fit.glmpca(y=Y, x=X, ncomp=ncomp, family=family, method = "fisher", verbose = FALSE)
  saveRDS(object = fit, file = paste0("Output_files/glmpca_200k_fisher_def_", i, ".RDS"))
  rm(fit)

}

#' 
