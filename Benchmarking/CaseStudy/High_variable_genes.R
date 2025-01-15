#' ---
#' title: "Untitled"
#' author: "Alex"
#' date: "2024-03-13"
#' output: html_document
#' ---
## -------------------------------------------------------------------------------------------------
library(HDF5Array)
library(scater)
library(scran)

tenx <- loadHDF5SummarizedExperiment(
                             dir =  "TENxBrainData_preprocessed_default", 
                             prefix="")

tenx <- logNormCounts(tenx)


#' 
## -------------------------------------------------------------------------------------------------
dec.sce <- modelGeneVar(tenx)

saveRDS(dec.sce, file = "Output_files/dec.sce.RDS")

#' 
