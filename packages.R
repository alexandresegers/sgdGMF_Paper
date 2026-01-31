
# Check whether a package is already installed, if not than install it
check.install.packages = function(pkg, type = c("default", "github", "bioc"), path = "", loadpkg = FALSE, ...) {
  type = match.arg(type)

  # Check if package is installed
  check = pkg %in% installed.packages()[,"Package"]
  cat("Already installed:", check, "\n")

  # If it is not, install it
  if (!check) {
    repo = if (path != "") paste(path, pkg, sep = "/") else pkg
    tryCatch({
      switch(type,
             "default" = utils::install.packages(pkg, ...),
             "github" = devtools::install_github(repo, ...),
             "bioc" = BiocManager::install(pkg, ...))
    }, error = function(e) {
      message("Failed to install ", pkg, ": ", e$message)
    })
  }

  # Load the package if not already loaded
  if (loadpkg) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      warning("Package ", pkg, " is installed but failed to load")
    }
  }
}

# Package managers
check.install.packages("devtools")
check.install.packages("BiocManager")

# sgdGMF
check.install.packages("sgdGMF")

# Fast linear algebra
check.install.packages("Matrix")
check.install.packages("irlba")
check.install.packages("svd")
check.install.packages("RSpectra")

# Generalized factor models
check.install.packages("gllvm")
check.install.packages("GFM")
check.install.packages("COAP")
check.install.packages("rrpack")

# Generalized PCA
check.install.packages("glmpca")
check.install.packages("PoissonPCA")
check.install.packages("gmf", type = "github", path = "kidzik")

# Matrix factorization
check.install.packages("NMF")
check.install.packages("cmfrec")
check.install.packages("NNLM", type = "github", path = "linxihui")

# Low-dimensional embedding
check.install.packages("Rtsne")
check.install.packages("umap")

# Omics classes and methods
check.install.packages("scry", type = "bioc")
check.install.packages("splatter", type = "bioc")
check.install.packages("scater", type = "bioc")
check.install.packages("zinbwave", type = "bioc")
check.install.packages("NewWave", type = "github", path = "fedeago")
# check.install.packages("NewWave", type = "bioc")

# Clustering
check.install.packages("cluster")
check.install.packages("bluster", type = "bioc")

# Visualization
check.install.packages("dplyr")
check.install.packages("ggplot2")
check.install.packages("reshape2")
check.install.packages("ggpubr")
check.install.packages("GGally")
check.install.packages("factoextra")

