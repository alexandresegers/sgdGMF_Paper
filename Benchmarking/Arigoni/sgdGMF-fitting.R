#' ---
#' title: "Model selection and hvg influence"
#' output: html_document
#' date: '2024-03-10'
#' ---
#' 
## --------------------------------------------------------------------------------------------------
source("utilities.R")
library(sgdGMF)
library(dplyr)

sim_full <- readRDS(file = "Data/BE1/final/sce.RDS")


#' 
## --------------------------------------------------------------------------------------------------
set.seed(100)

## LOAD DATA ----


hvg <- (rownames(sim_full)[order(rowData(sim_full)$biological.var, decreasing = TRUE)])[1:500]
sim <- sim_full[hvg,]
logcounts = as.data.frame(as.matrix(logcounts(sim)))
counts = as.data.frame(as.matrix(counts(sim)))
cells = as.data.frame(colData(sim))
genes = as.data.frame(rowData(sim))
meta = metadata(sim)
groups = as.numeric(as.factor(cells$celltype))

n = ncol(sim)
m = nrow(sim)

X = model.matrix(~ 1, data = cells)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)


family = poisson()
method = "sgd"
sampling = "sampling"
control = list("parallel" = TRUE, "nthreads" = 8)

gmf.fit.list = list()
ncomp <- c(1:10, 15, 20, 30, 40, 50)
control.init = list(method = "ols", type = "link")
control.alg = list(maxiter = 1000, size = c(100,100), frequency = 250)

for (i in 1:length(ncomp)) {
  cat(" Rank =", ncomp[i], "\n")
  # Fit the model
  fit = sgdGMF::sgdgmf.fit(Y, X, Z, ncomp = ncomp[i], family = family,
                           method = method, sampling = sampling, 
                           control.init = control.init, 
                           control.alg = control.alg)
  saveRDS(object = fit, file = paste0("Output_models/sgdGMF-B-SGD-fit-", ncomp[i], ".RDS"))
  
  # Store the estimated model
  gmf.fit.list[[i]] = fit
}

saveRDS(object = gmf.fit.list, file = paste0("Output_models/sgdGMF-B-SGD-fit-list.RDS"))
  
  

#' 
#' 
