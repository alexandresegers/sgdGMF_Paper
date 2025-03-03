#' ---
#' title: "Untitled"
#' output: html_document
#' date: '2024-05-21'
#' ---
#' 
## ----setup, include=FALSE--------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## --------------------------------------------------------------------------------------------------
library(bluster)
library(SummarizedExperiment)
library(NewWave)
library(scater)
library(pheatmap)
library(igraph)
library(scran)

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

source("utilities.R")

set.seed(100)

#' 
## --------------------------------------------------------------------------------------------------
sce <- readRDS(file = "Data/BE1/final/sce.RDS")
sce <- sce[(rownames(sce)[order(rowData(sce)$biological.var, 
                                decreasing = TRUE)])[1:500],]
celltypes <- colData(sce)$celltype

saveRDS(celltypes, file = "Output_models/celltypes.RDS")

sgdmodels_cv <- readRDS(file = "Output_models/sgdGMF-B-SGD-cv.RDS")
sgd_models <- readRDS(file = "Output_models/sgdGMF-B-SGD-fit-list.RDS")
eigenvalues <- readRDS(file = "Output_models/eigenvalues.RDS")



#' 
#' # Model selection
#' 
#' ## Selection criteria
## --------------------------------------------------------------------------------------------------
df_model_selection <- data.frame(vars = rep(c("AIC","BIC","Deviance"), 
                                            each = nrow(sgdmodels_cv$summary.cv)),
                 values = c(sgdmodels_cv$summary.cv$aic, sgdmodels_cv$summary.cv$bic,
                            sgdmodels_cv$summary.cv$dev),
                 ncomp = as.factor(rep(sgdmodels_cv$summary.cv$ncomp, times = 3)))

df_eigenvalues <- data.frame("eigenvalues" = eigenvalues$lambdas, 
                             "ncomp" = c(1:500), facet = "1")

saveRDS(object = df_model_selection, file = "Output_models/df_model_selection.RDS")
saveRDS(object = df_eigenvalues, file = "Output_models/df_eigenvalues.RDS")

#' 
#' ## Purity
## --------------------------------------------------------------------------------------------------
dimensions <- c(1:10, 15, 20, 30, 40, 50)
purity_rescaled <- lapply(1:15, function(x) {
  U_rescaled <- prcomp(sgd_models[[x]]$U %*% t(sgd_models[[x]]$V), rank = dimensions[x])$x
  pur <- neighborPurity(U_rescaled, clusters = celltypes)
  pur$dim <- c(dimensions[x])
  return(pur)
  })

purity_rescaled_merged <- do.call(rbind, purity_rescaled)
purity_rescaled_merged$dim <- as.factor(purity_rescaled_merged$dim)
purity_rescaled_mean <- as.data.frame(purity_rescaled_merged) %>% group_by(maximum, dim) %>% 
  summarise(mean = mean(purity))
purity_rescaled_mean$facet <- "1"

saveRDS(object = purity_rescaled_mean, file = "Output_models/purity_rescaled_mean.RDS")


#' 
#' ## t-sne plots
#' 
## --------------------------------------------------------------------------------------------------
tsne_models <- lapply(c(1:15), function(x){
  if(x == 1){
    tsne <- scale(Rtsne::Rtsne(as.matrix(sgd_models[[1]]$U, ncol = 1), dims = 2, 
                                 verbose = TRUE, num_threads = 10)$Y)
  return(tsne)
  }
  tsne <- scale(Rtsne::Rtsne(sgd_models[[x]]$U, dims = 2, 
                                 verbose = TRUE, num_threads = 10)$Y)
  return(tsne)
})

df_tsne_plot <- as.data.frame(do.call(rbind,tsne_models))
df_tsne_plot$dim <- rep(c(1:10, 15, 20, 30, 40, 50), each = nrow(tsne_models[[1]]))
df_tsne_plot$celltypes <- rep(celltypes, times = 15)

saveRDS(object = df_tsne_plot, file = "Output_models/df_tsne_plot.RDS")



#' 
#' 
#' ## Clustering t-sne
#' 
## --------------------------------------------------------------------------------------------------
set.seed(100)
tsne_9 <- tsne_models[[9]]
colnames(tsne_9) <- c("T1","T2")
scores <- prcomp(sgd_models[[9]]$U %*% t(sgd_models[[9]]$V),rank = 9)$x
reducedDim(sce, "SGD_rescaled_9") <- scores
graph_9 <- buildSNNGraph(x = sce, 
                       use.dimred = "SGD_rescaled_9") 

resolution_param <- c(0.0038)

leiden_clustering_9 <- igraph::cluster_leiden(graph_9, resolution_parameter = resolution_param)

set.seed(100)
tsne <- tsne_models[[11]]
colnames(tsne) <- c("T1","T2")
scores <- prcomp(sgd_models[[11]]$U %*% t(sgd_models[[11]]$V),rank = 15)$x
reducedDim(sce, "SGD_rescaled_15") <- scores
graph <- buildSNNGraph(x = sce, 
                       use.dimred = "SGD_rescaled_15") 


resolution_param <- c(0.004)

leiden_clustering_15 <- igraph::cluster_leiden(graph, resolution_parameter = resolution_param)


set.seed(100)
tsne_30 <- tsne_models[[13]]
colnames(tsne_30) <- c("T1","T2")
scores <- prcomp(sgd_models[[13]]$U %*% t(sgd_models[[13]]$V),rank = 30)$x
reducedDim(sce, "SGD_rescaled_30") <- scores
graph_30 <- buildSNNGraph(x = sce, 
                       use.dimred = "SGD_rescaled_30") 


resolution_param <- c(0.004)

leiden_clustering_30 <- igraph::cluster_leiden(graph_30, resolution_parameter = resolution_param)


#' 
## --------------------------------------------------------------------------------------------------
df_clust_tsne <- data.frame("T1" = c(tsne_9[,1], tsne[,1], tsne_30[,1]), 
                            "T2" = c(tsne_9[,2], tsne[,2], tsne_30[,2]),
                            "membership" = c(factor(leiden_clustering_9$membership, 
                                                   levels = c(1, 2, 4, 6, 5, 3, 7),
                                                    labels = c("1","2","3","4","5","6","7")),
                                             factor(leiden_clustering_15$membership,
                                                    levels = c(1, 2, 6, 5, 4, 3, 7),
                                                    labels = c("1","2","3","4","5","6","7")),
                                             factor(leiden_clustering_30$membership,
                                                    levels = c(1, 5, 6, 4, 3, 2, 7),
                                                    labels = c("1","2","3","4","5","6","7"))),
                            "dims" = as.factor(rep(c(9, 15, 30), each = 26395)))

saveRDS(object = df_clust_tsne, file = "Output_models/df_clust_tsne.RDS")


#' 
#' ## clustering confusion matrix
## --------------------------------------------------------------------------------------------------
tab_9 <- table(sce$celltype, factor(leiden_clustering_9$membership, 
                                                  levels = c(1, 2, 4, 6, 5, 3, 7),
                                                    labels = c("1","2","3","4","5","6","7")))

tab_15 <- table(sce$celltype, factor(leiden_clustering_15$membership,
                                                    levels = c(1, 2, 6, 5, 4, 3, 7),
                                                    labels = c("1","2","3","4","5","6","7")))

tab_30 <- table(sce$celltype, factor(leiden_clustering_30$membership,
                                                    levels = c(1, 5, 6, 4, 3, 2, 7),
                                                    labels = c("1","2","3","4","5","6","7")))



#' 
## --------------------------------------------------------------------------------------------------
df_tile_confusion <- data.frame("hits" = c(tab_9, tab_15, tab_30), 
                           "celltypes" = factor(rownames(tab_9), 
                                                levels = c("A549", "CCL-185-IG",
                                                           "CRL5868", "DV90",
                                                           "HCC78", "HTB178", "PC9")),
                           "cluster" = rep(c(1:7), each = 7),
                           "dims" = factor(rep(c(9, 15, 30), each = 49)))

saveRDS(object = df_tile_confusion, file = "Output_models/df_tile_confusion.RDS")


#' 
#' 
#' # Supplementary Figure: different HVG
#' 
## --------------------------------------------------------------------------------------------------
# set.seed(100)
# dimensions <- c(1:10, 15, 20, 30, 40, 50)
# 
# df_purity_hvg <- data.frame("celltype" = NULL,"dim" = NULL, "mean" = NULL, "hvg" = NULL)
# for (i in c(100,200,500,1000,1500,2000)){
#   if(i == 500){
#     sgd_model <- readRDS(file = paste0("Output_models/sgdGMF-B-SGD-fit-list.RDS"))
# 
#   }else{
#   sgd_model <- readRDS(file = paste0("Output_models/sgdGMF-B-SGD-fit-list-",i,"-hvg.RDS"))
#   }
#   
# 
# 
#     purity_rescaled_hvg <- lapply(1:15, function(x) {
#       U_rescaled <- prcomp(sgd_model[[x]]$U %*% t(sgd_model[[x]]$V))$x
#       pur <- neighborPurity(U_rescaled, clusters = celltypes)
#       pur$dim <- c(dimensions[x])
#       return(pur)
#       })
#     purity_rescaled_merged <- do.call(rbind, purity_rescaled_hvg)
#     purity_rescaled_merged$dim <- as.factor(purity_rescaled_merged$dim)
#     purity_rescaled_mean <- as.data.frame(purity_rescaled_merged) %>% group_by(maximum, dim) %>%
#       summarise(mean = mean(purity))
#     
#     df_purity_hvg <- rbind(df_purity_hvg, data.frame(purity_rescaled_mean$maximum, 
#                               purity_rescaled_mean$dim,
#                               purity_rescaled_mean$mean, 
#                               rep(i, nrow(purity_rescaled_mean))))
#   
# }
# colnames(df_purity_hvg) <- c("celltype", "dim", "mean", "hvg")
# df_purity_hvg$dim <- as.factor(df_purity_hvg$dim)
# df_purity_hvg$hvg <- as.factor(df_purity_hvg$hvg)
# df_purity_hvg$celltype <- as.factor(df_purity_hvg$celltype)
# 
# saveRDS(object = df_purity_hvg, file = "Output_models/df_purity_hvg.RDS")


#' 
#' 
## --------------------------------------------------------------------------------------------------
# set.seed(100)
# df_tsne_hvg <- data.frame("T1" = NULL,"T2" = NULL, "dim" = NULL, "hvg" = NULL, "celltypes" = NULL)
# for (i in c(100,200,500,1000,1500,2000)){
#   if(i == 500){
#     sgd_model <- readRDS(file = paste0("Output_models/sgdGMF-B-SGD-fit-list.RDS"))
# 
#   }else{
#   sgd_model <- readRDS(file = paste0("Output_models/sgdGMF-B-SGD-fit-list-",i,"-hvg.RDS"))
#   }
#   
#   for(j in c(5,10,13)){
#     set.seed(100)
#     tsne <- scale(Rtsne::Rtsne(sgd_model[[j]]$U, dims = 2, 
#                                  verbose = TRUE, num_threads = 8)$Y)
#     df_tsne_hvg <- rbind(df_tsne_hvg, data.frame(tsne[,1], 
#                               tsne[,2],
#                               rep(j, nrow(sgd_model[[j]]$U)), 
#                               rep(i, nrow(sgd_model[[j]]$U)), 
#                               celltypes))
#   }
# }
# colnames(df_tsne_hvg) <- c("T1", "T2", "dim", "hvg", "celltypes")
# df_tsne_hvg$dim[df_tsne_hvg$dim == 13] <- 30
# df_tsne_hvg$dim <- as.factor(df_tsne_hvg$dim)
# df_tsne_hvg$hvg <- as.factor(df_tsne_hvg$hvg)
# df_tsne_hvg$celltypes <- as.factor(df_tsne_hvg$celltypes)
# 
# saveRDS(object = df_tsne_hvg, file = "Output_models/df_tsne_hvg.RDS")


#' 
#' 
#' # Supplementary fig: NewWave comparison:
#' 
## --------------------------------------------------------------------------------------------------
matrix.deviance <- function (mu, y, family = poisson())
{
dev = pointwise.deviance(mu, y, family)
dev = sum(dev, na.rm = TRUE)
return(dev)
}
pointwise.deviance <- function(mu, y, family = poisson()){
if (length(mu) == 1) {
mut = y
mut[] = mu
mu = mut
}
  nona = !is.na(y)
dev = y
dev[] = NA
dev[nona] = family$dev.resids(y[nona], mu[nona], 1)
return(dev)
}

#' 
#' 
## --------------------------------------------------------------------------------------------------
models <- readRDS(file = "Output_models/NewWave_comparison_results_train_test_15_ncomp.RDS")
train_test <- readRDS(file = "Output_models/train_test_15_ncomp.RDS")

model.nbwave <- models$NewWave
model.bsgd_poisson <- models$`Poisson-BSGD`
model.bsgd <- models$`NB-BSGD`
model.avagrad <- models$Avagrad
model.fisher <- models$Fisher


test = train_test$test
train = train_test$train

models_full <- readRDS(file = "Output_models/NewWave_comparison_results_all_data_15_ncomp.RDS")


modelfull.nbwave <- models_full$NewWave
modelfull.bsgd_poisson <- models_full$`Poisson-BSGD`
modelfull.bsgd <- models_full$`NB-BSGD`
modelfull.avagrad <- models_full$Avagrad
modelfull.fisher <- models_full$Fisher

#' 
## --------------------------------------------------------------------------------------------------
m <- nrow(sce)
n <- ncol(sce)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts(sce))

set.seed(100)
full.sample.error = error.summary(
  Y, 
  as.factor(celltypes),
  model.nbwave, 
  model.bsgd,
  model.bsgd_poisson,
  model.avagrad,
  model.fisher)
print(full.sample.error)

# Out-of-sample error
in.sample.error = error.summary(
  train, 
  as.factor(celltypes),
  model.nbwave, 
  model.bsgd,
  model.bsgd_poisson,
  model.avagrad,
  model.fisher)
print(in.sample.error)

# Out-of-sample error
out.sample.error = error.summary(
  test, 
  as.factor(celltypes),
  model.nbwave,
  model.bsgd,
  model.bsgd_poisson,
  model.avagrad,
  model.fisher)
print(out.sample.error)

df_errors = data.frame(rbind(full.sample.error, in.sample.error, out.sample.error))
df_errors$Sample = rep(c("full", "train", "test"), each = 5)
df_errors$Model = factor(rep(c("NewWave","SGD-NB","SGD-Poisson", "Avagrad", "Fisher"),times = 3),
                  levels= c("Avagrad", "Fisher", "NewWave", "SGD-NB", "SGD-Poisson"))
df_errors$Sample = as.factor(df_errors$Sample)
df_errors$Time = as.numeric(df_errors$Time)
df_errors$RSS = as.numeric(df_errors$RSS)
df_errors$Cos = as.numeric(df_errors$Cos)
df_errors$Dev = as.numeric(df_errors$Dev)

saveRDS(df_errors, file = "Output_models/df_errors.RDS")

#' 
## --------------------------------------------------------------------------------------------------
set.seed(100)
purity_NewWave <- neighborPurity(modelfull.nbwave$u%*%diag(modelfull.nbwave$d), 
                                 clusters = celltypes)
purity_sgd_nb <- neighborPurity(modelfull.bsgd$u%*%diag(modelfull.bsgd$d), clusters = celltypes)

purity_sgd_poisson <- neighborPurity(modelfull.bsgd_poisson$u%*%diag(modelfull.bsgd_poisson$d),
                                     clusters = celltypes)
purity_avagrad <- neighborPurity(modelfull.avagrad$u%*%diag(modelfull.avagrad$d), clusters = celltypes)

purity_fisher <- neighborPurity(modelfull.fisher$u%*%diag(modelfull.fisher$d), clusters = celltypes)

purity_newwave_comparison <- rbind(purity_NewWave, purity_sgd_nb, purity_sgd_poisson, 
                                   purity_avagrad, purity_fisher)
purity_newwave_comparison$method <- factor(rep(c("NewWave","SGD-NB","SGD-Poisson", "Avagrad", "Fisher"),
                                                  each = 26395), levels= c("Avagrad", "Fisher", "NewWave", "SGD-NB", "SGD-Poisson"))


purity_newwave_comparison$maximum <- as.factor(purity_newwave_comparison$maximum)

purity_newwave_comparison_mean <- as.data.frame(purity_newwave_comparison) %>% 
  group_by(maximum, method) %>%
  summarise(mean = mean(purity))

saveRDS(purity_newwave_comparison_mean, file = "Output_models/purity_newwave_comparison_mean.RDS")


#' 
#' 
## --------------------------------------------------------------------------------------------------
df_tsne_newwave_comparison <- data.frame("T1" = c(modelfull.nbwave$tsne[,1], modelfull.bsgd$tsne[,1],
                                                modelfull.bsgd_poisson$tsne[,1], modelfull.avagrad$tsne[,1],
                                                  modelfull.fisher$tsne[,1]),
                                         "T2" = c(modelfull.nbwave$tsne[,2], modelfull.bsgd$tsne[,2],
                                                modelfull.bsgd_poisson$tsne[,2], modelfull.avagrad$tsne[,2],
                                                  modelfull.fisher$tsne[,2]),
                                         "method" = factor(rep(c("NewWave", "SGD-NB", "SGD-Poisson",
                                                            "Avagrad", "Fisher"), each = 26395),
                                                           levels= c("Avagrad", "Fisher", "NewWave", 
                                                                     "SGD-NB", "SGD-Poisson")),
                                         "celltype" = rep(celltypes, times = 5))


saveRDS(df_tsne_newwave_comparison, file = "Output_models/df_tsne_newwave_comparison.RDS")


#' 
