#' ---
#' title: "Untitled"
#' author: "Alex"
#' date: "2024-03-06"
#' output: html_document
#' ---
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(scDblFinder)


#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
counts <- readRDS(file = "Data/BE1/filtered/counts.RDS")
barcodes <- readRDS(file = "Data/BE1/filtered/barcodes.RDS")
celltype <- readRDS( file = "Data/BE1/filtered/celltype.RDS")
features <- readRDS(file = "Data/BE1/filtered/features.RDS")

#' 
#' 
#' 
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(100)
colnames(counts) <- barcodes
rownames(counts) <- features[,1]

colnames(features) <- c("ensembl_name", "gene_name","type_of_expression")


sce <- SingleCellExperiment::SingleCellExperiment(assays = list("counts" = counts),
                                                  rowData = features,
                                                  colData = DataFrame("celltype" = celltype))

sce <- sce[features[,3] == "Gene Expression",]
sce <- sce[,celltype != "PBMCs"]

#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------

is.mito <- grep(pattern = "^MT-", x = features[,2])

library(scuttle)
df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))

colData(sce) <- cbind(colData(sce), df)

#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
reasons <- perCellQCFilters(df, sub.fields="subsets_Mito_percent")
colSums(as.matrix(reasons))

#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
summary(reasons$discard)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------
sce$discard <- reasons$discard


#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(data = as.data.frame(colData(sce)), aes( x= "", y = subsets_Mito_percent, col = discard)) +
    geom_jitter()

ggplot(data = as.data.frame(colData(sce)), aes( x= "", y = total, col = discard)) +
    geom_jitter()

ggplot(data = as.data.frame(colData(sce)), aes( x= "", y = detected, col = discard)) +
    geom_jitter()

#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
sce <- sce[,!sce$discard]
sce


#' 
#' # Features
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
sce <- logNormCounts(sce)
sce

#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
dec.sce <- modelGeneVar(sce)
fit.sce <- metadata(dec.sce)

plot(fit.sce$mean, fit.sce$var, xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curve(fit.sce$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

rowData(sce)$biological.var <- dec.sce$bio


#' 
#' 
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
hvg.sce.var <- getTopHVGs(dec.sce, n=500)
head(hvg.sce.var)
sce <- runPCA(sce, subset_row = hvg.sce.var)


#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(100)

dbl.dens <- computeDoubletDensity(sce, subset.row=hvg.sce.var,
                                  d=ncol(reducedDim(sce)))
summary(dbl.dens)

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call")
summary(dbl.calls)

#' 
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
sce <- sce[,!(dbl.calls == "doublet")]

#' 
#' 
#' 
#' # Rerunning some characteristics on final data
#' 
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
sce <- logNormCounts(sce)
dec.sce <- modelGeneVar(sce)
fit.sce <- metadata(dec.sce)

plot(fit.sce$mean, fit.sce$var, xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curve(fit.sce$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

rowData(sce)$biological.var <- dec.sce$bio

#' 
#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
hvg.sce.var <- getTopHVGs(dec.sce, n=500)
head(hvg.sce.var)
sce <- runPCA(sce, subset_row = hvg.sce.var)


#' 
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(object = sce,
        file = "Data/BE1/final/sce.RDS")
save.image(file = "Data/BE1/final/sce.RData")

