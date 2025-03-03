#' ---
#' title: "Untitled"
#' output: html_document
#' date: '2024-03-15'
#' ---
#' 
## ----setup, include=FALSE--------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' 
## --------------------------------------------------------------------------------------------------
library(sgdGMF)
library(dplyr)

#' 
## --------------------------------------------------------------------------------------------------
set.seed(100)

sim_full <- readRDS(file = "Data/BE1/final/sce.RDS")

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
sampling = "block"
control.init = list(method = "ols", type = "link")
control.alg = list(maxiter = 1000, size = c(100,25), frequency = 250)
control.cv = list(nfolds = 5)

gmf.fit = sgdgmf.cv(Y, X, Z, ncomps = c(1:10,15,20,30,40,50), family = family,
                    method = method, sampling = sampling,
                    control.init = control.init, control.alg = control.alg, control.cv = control.cv)

saveRDS(gmf.fit, file = "Output_models/sgdGMF-B-SGD-cv.RDS")

