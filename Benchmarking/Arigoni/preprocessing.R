
# Clean the workspace
rm(list = ls())
graphics.off()

# Libraries
library(Matrix)

# Global variables
SHOW = TRUE
SAVE = TRUE
THR = 0.1

# Load data function
load.data = function () {
  
  thr = THR
  
  data.list = list()
  cell.list = c("A549", "CCL-185-IG", "CRL5868", "DV90", "HCC78", "HTB178", "PBMCs", "PC9")
  
  H = length(cell.list)
  for (h in 1:H) {
      
    folder = cell.list[h]
    barcodes = paste("Data/BE1/raw", folder, "barcodes.tsv.gz", sep = "/")
    features = paste("Data/BE1/raw", folder, "features.tsv.gz", sep = "/")
    metrics = paste("Data/BE1/raw", folder, "metrics_summary.csv", sep = "/")
    matrix = paste("Data/BE1/raw", folder, "matrix.mtx.gz", sep = "/")
    
    cat(" - cell:", folder, "\n")
    data = list()
    data$cell = folder
    data$barcodes = as.matrix(readr::read_tsv(gzfile(barcodes), col_names = FALSE, show_col_types = FALSE))
    data$features = as.matrix(readr::read_tsv(gzfile(features), col_names = FALSE, show_col_types = FALSE))
    data$metrics = read.csv(metrics, header = TRUE)
    data$counts = Matrix::readMM(gzfile(matrix))
    data$gene.means = rowMeans(data$counts)
    data$cell.means = colMeans(data$counts)
    data$gene.stds = apply(data$counts, 1, sd)
    data$cell.stds = apply(data$counts, 2, sd)
    data.list[[h]] = data
    rm(data); gc()
  }
  
  gene.means = rowMeans(do.call(cbind, lapply(data.list, function (x) x$gene.means)))
  gene.stds = rowMeans(do.call(cbind, lapply(data.list, function (x) x$gene.stds)))
  
  idx = which(gene.means / mean(gene.means) > thr)
  
  cat(" - gene: selection \n")
  for (h in 1:H) {
    data = data.list[[h]]
    data$celltype = rep(data$cell, times = ncol(data$counts))
    data$features = data$features[idx,]
    data$counts = data$counts[idx,]
    data$gene.means = data$gene.means[idx]
    data$gene.stds = data$gene.stds[idx]
    data.list[[h]] = data
    rm(data); gc()
  }
  
  return (data.list)
}

# Load the data
data.list = load.data()

# separate the lists of matrices, features, barcodes and cell-types
mat.list = lapply(data.list, function (x) x$counts)
ftr.list = lapply(data.list, function (x) x$features)
bar.list = lapply(data.list, function (x) x$barcodes)
cll.list = lapply(data.list, function (x) x$celltype)

# Concatenate the information to get the final data matrices
mat = do.call(cbind, mat.list)
ftr = do.call(cbind, ftr.list)
bar = do.call(c, bar.list)
cll = do.call(c, cll.list)

if (SAVE) {
  # Save the preprocessed data
  writeMM(mat, file = "Data/BE1/filtered/counts.mtx")
  write.table(bar, file = "Data/BE1/filtered/barcodes.csv", sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(cll, file = "Data/BE1/filtered/celltype.csv", sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
  write.table(ftr[,1:3], file = "Data/BE1/filtered/features.csv", sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
}

SAVERDS = TRUE
if(SAVERDS){
    saveRDS(mat, file = "Data/BE1/filtered/counts.RDS")
    saveRDS(bar, file = "Data/BE1/filtered/barcodes.RDS")
    saveRDS(cll, file = "Data/BE1/filtered/celltype.RDS")
    saveRDS(ftr[,1:3], file = "Data/BE1/filtered/features.RDS")
    
}

if (SHOW) {
  # Sort the genes by their mean
  idx = order(rowMeans(mat), decreasing = TRUE)
  mat = t(mat[idx, ])
  
  # Project the cells into a 2D space via PCA and tSNE
  X = log1p(scale(mat[, 1:1000], center = FALSE, scale = TRUE))
  pca = RSpectra::svds(X, 2)
  tsn = Rtsne::Rtsne(X, 2, num_threads = 8, verbose = TRUE)
  
  plot(pca$u, col = as.numeric(as.factor(cll)))
  plot(tsn$Y, col = as.numeric(as.factor(cll)))
}


# counts = Matrix::readMM("counts.mtx")
# genetypes = c(read.table(file = "barcodes.csv", sep = ",", dec = ".", col.names = FALSE))
# celltypes = c(read.table(file = "celltype.csv", sep = ",", dec = ".", col.names = FALSE))
# save(counts, genetypes, celltypes, file = "BE1.RData")
