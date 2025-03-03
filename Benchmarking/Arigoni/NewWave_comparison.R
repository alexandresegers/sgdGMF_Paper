## --------------------------------------------------------------------------------------------------
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load the package
#devtools::load_all()

## Load the utility functions
source("utilities.R")
library(MASS)
library(sgdGMF)
library(dplyr)

#' 
#' 
## --------------------------------------------------------------------------------------------------
set.seed(100)
## DATA SIMULATION ----
SAVE = FALSE
SHOW = TRUE


sim_full <- readRDS(file = "Data/BE1/final/sce.RDS")
hvg <- (rownames(sim_full)[order(rowData(sim_full)$biological.var,
                                   decreasing = TRUE)])[1:500]
sim <- sim_full[hvg,]

assays(sim)$counts <- as.matrix(assays(sim)$counts)

m <- nrow(sim)
n <- ncol(sim)

## DATA EXTRACTION ----
logcounts = as.data.frame(as.matrix(logcounts(sim)))
counts = as.data.frame(as.matrix(counts(sim)))
cells = as.data.frame(colData(sim))
genes = as.data.frame(rowData(sim))
meta = metadata(sim)


## TRAIN-TEST SPLIT ----
X = model.matrix(~ 1, data = cells)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)

ncomp = 15
family = poisson()

# Sparsification and matrix completion
data = train.test.split(Y, test = 0.3)

test = Y
train = Y
ctrain = Y

test[data$test] = NA
train[data$train] = NA
ctrain = naive.completion(train)

saveRDS(list("train" = train, "test" = test), file = "Output_models/train_test_15_ncomp.RDS")



model.bsgd_poisson = fit.C.bsgd(y = train, x = X, z = Z, ncomp = ncomp, 
                                family = family, verbose = FALSE, maxiter = 1000)
model.bsgd = fit.C.bsgd(y = train, x = X, z = Z, ncomp = ncomp, family = neg.bin(1), 
                        verbose = FALSE, maxiter = 1000)
model.nbwave = fit.nbwave(y = ctrain, x = X, z = Z, ncomp = ncomp, family = neg.bin(1), 
                          verbose = FALSE, maxiter = 100, tol = 1e-04)
avagrad <- fit.glmpca(y=ctrain, x=X, z=Z, ncomp=ncomp, family=family, method = "avagrad", verbose= FALSE, maxiter=1000, tol=1e-04)
fisher <- fit.glmpca(y=ctrain, x=X, z=Z, ncomp=ncomp, family=family, method = "fisher", verbose= FALSE, maxiter=200, tol=1e-05)

## MODEL CHECK ----

# Execution time
{
cat("NBWaVE: ", model.nbwave$time[3], "s \n")
cat("B-SGD:  ", model.bsgd$time[3], "s \n")
cat("B-SGD-Poisson:  ", model.bsgd_poisson$time[3], "s \n")
}

results <- list("NewWave" = model.nbwave, "NB-BSGD" = model.bsgd, 
                "Poisson-BSGD" = model.bsgd_poisson, 
                "Avagrad" = avagrad, "Fisher" = fisher)
saveRDS(results, file = "Output_models/NewWave_comparison_results_train_test_15_ncomp.RDS")




model.bsgd_poisson = fit.C.bsgd(y = Y, x = X, z = Z, ncomp = ncomp, family = family, 
                                verbose = FALSE, maxiter = 1000)
model.bsgd = fit.C.bsgd(y = Y, x = X, z = Z, ncomp = ncomp, family = neg.bin(1), 
                        verbose = FALSE, maxiter = 1000)
model.nbwave = fit.nbwave(y = Y, x = X, z = Z, ncomp = ncomp, family = neg.bin(1), 
                          verbose = FALSE, maxiter = 100, tol = 1e-04)
avagrad <- fit.glmpca(y=Y, x=X, z=Z, ncomp=ncomp, family=family, method = "avagrad", 
                      verbose= FALSE, maxiter=1000, tol=1e-04)
fisher <- fit.glmpca(y=Y, x=X, z=Z, ncomp=ncomp, family=family, method = "fisher",
                     verbose= FALSE, maxiter=200, tol=1e-05)

## MODEL CHECK ----

# Execution time
{
cat("NBWaVE: ", model.nbwave$time[3], "s \n")
cat("B-SGD:  ", model.bsgd$time[3], "s \n")
cat("B-SGD-Poisson:  ", model.bsgd_poisson$time[3], "s \n")
}

results <- list("NewWave" = model.nbwave, "NB-BSGD" = model.bsgd, 
                "Poisson-BSGD" = model.bsgd_poisson, "Avagrad" = avagrad, 
                "Fisher" = fisher)
saveRDS(results, file = "Output_models/NewWave_comparison_results_all_data_15_ncomp.RDS")


#' 
#' 
