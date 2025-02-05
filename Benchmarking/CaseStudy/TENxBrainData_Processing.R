#' 
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, 
                      error = FALSE, message = FALSE, 
                      cache = FALSE)

#' 
## ----create-folders-------------------------------------------------------------------------------
#Create the necessary folder structure

library(doParallel)
ncores <- 5
doParallel::registerDoParallel(ncores)

#' 
#' ## Data import 
## ----data-import----------------------------------------------------------------------------------
library(TENxBrainData)
tenx <- TENxBrainData()
tenx

## Check that the counts object is a HDF5Array
seed(counts(tenx))
assay(tenx)

#' 
#' Run the next chunk to pick a subset of cells to debug code.
#' 
## ----tmp------------------------------------------------------------------------------------------
# set.seed(1234)
# keep_cells <- sample(1:ncol(tenx),size = 1000)
# tenx <- tenx[,keep_cells]

#' 
#' # Quality control and 
#' 
#' ## Removing low-quality cells
#' 
#' First, we use the `scater` package to compute a set of 
#' QC measures and filter out the low-quality samples.
#' 
## ----qc-calculate---------------------------------------------------------------------------------
library(scater)
system.time(tenx <- addPerCellQC(tenx, 
                          subsets = list(Mito = grep("^mt", rowData(tenx)$Symbol)),
                          BPPARAM = BiocParallel::MulticoreParam(ncores)))

#' 
#' We remove cells with high proportion of mitocondrial 
#' reads, using it as a proxy for cell damage. 
#' 
## ----qc-cell-filter-------------------------------------------------------------------------------
high_mito <- isOutlier(tenx$subsets_Mito_percent, 
                       nmads = 3, type="higher")
table(high_mito)

tenx <- tenx[,!high_mito]
tenx

#' 
#' ## Removing lowly expressed genes
#' 
#' Next, we remove the lowly expressed genes. Here, 
#' we keep only those genes that have at least 1 UMI 
#' in at least 1% of the data. These threshold are
#' dataset-specific and may need to be taylored to 
#' specific applications.
#' 
## ----qc-gene-filter-------------------------------------------------------------------------------
num_reads <- 1
num_cells <- 0.01*ncol(tenx)
system.time(keep <- which(DelayedArray::rowSums(counts(tenx) >= num_reads ) >= num_cells))
tenx <- tenx[keep,]
tenx

#' 
#' These leaves us with `length(keep)` genes.
#' 
#' 
#' Save tenx object to
## ----save-hdf5-files------------------------------------------------------------------------------
library(HDF5Array)

saveHDF5SummarizedExperiment(tenx,
                             dir =  "TENxBrainData_preprocessed",
                             prefix="", replace=TRUE,
                             chunkdim=c(nrow(tenx),1),
                             level=NULL, verbose=FALSE)

saveHDF5SummarizedExperiment(tenx, 
                             dir =  "TENxBrainData_preprocessed_default", 
                             prefix="", replace=TRUE, 
                             level=NULL, verbose=FALSE)

