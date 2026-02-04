# ----------- Simulation benchmark

path_sim <- "Benchmarking/Simulation"

# First, all data is simulated with splatter, and all methods are run on these
# simulated datasets.
source(paste(path_sim, "simulation_MAIN.R", sep="/"))
source(paste(path_sim, "simulation_NB.R", sep="/"))

# Then, figures are made for the summary statistics and the tSNE projections.
source(paste(path_sim, "plot_sim_example_MAIN.R", sep="/"))

# Similarly, figures are made for the summary statistics in function of the
# matrix dimensions and latent space rank.
source(paste(path_sim, "plot_sim_summary_MAIN.R", sim="/"))
source(paste(path_sim, "plot_sim_summary_NB.R", sim="/"))


# ----------- Arigoni benchmark

path_arigoni <- "Benchmarking/Arigoni"

# First, the raw data of all celltypes is merged together in one file, and
# an initial quality control is done.
source(paste(path_arigoni, "preprocessing.R", sep="/"))

# Then, the data is additionally filtered. For example, PBMC cells are not used
# and cells with many mitochondrial reads are removed.
source(paste(path_arigoni, "Data_preparation.R", sep="/"))

# We calculate a scree-plot based on the eigenvalues for the model selection.
source(paste(path_arigoni, "Eigenvalues_model_selection.R", sep="/"))

# We also perform cross-validation based on model selection criteria and
# out of sample deviances to select the number of latent factors.
source(paste(path_arigoni, "sgdGMF-cv.R", sep="/"))

# sgdGMF estimation is performed using the optimal number of latent factors
# chosen in the cross-validation.
source(paste(path_arigoni, "sgdGMF-fitting.R", sep="/"))

# Benchmarking with glmPCA and NewWave is done. Here, these models are
# used to perform dimensionality reduction on the Arigoni dataset.
source(paste(path_arigoni, "NewWave_comparison.R", sep="/"))

# We here group all results in data.frames that can be used to make the plots
# of the manuscript.
source(paste(path_arigoni, "Figure_preparation.R", sep="/"))

# We here construct the plots of the manuscript.
source(paste(path_arigoni, "Figures_construction.R", sep="/"))



# ----------- CaseStudy benchmark

path_casestudy <- "Benchmarking/CaseStudy"

# First, the dataset is downloaded and quality control is performed to
# remove aberrant samples or genes.
source(paste(path_casestudy, "TENxBrainData_Processing.R", sep="/"))

# The variance of each gene is estimated to select the most variable
# genes of the dataset. These genes will be further used in the decomposition.
source(paste(path_casestudy, "High_variable_genes.R", sep="/"))

# We calculate a scree-plot based on the eigenvalues for the model selection.
source(paste(path_casestudy, "Model_selection_eigenvalues.R", sep="/"))

# sgdGMF estimation is performed using the optimal number of latent factors
# estimated with the scree plot.
source(paste(path_casestudy, "CaseStudy_Model_fitting_full.R", sep="/"))

# For the comparisons of the time- and memory-usage, sgdGMF is also run on 
# smaller parts of the data. Here, 5 times, the time is measured to compute the
# dimensionality reduction for 100.000, 200.000 and 300.000 cells in this
# dataset. This is also done for glmPCA with Avagrad estimation and with Fisher
# estimation. Also, NewWave is used on the subsampled dataset. Note that for
# reference, the COAP scripts are added, but are not run here due to errors
# that are returned.

path_casestudy_memory <- "Benchmarking/CaseStudy/Memory"

for (i in c(1:5)){
  
  script_name <- "CaseStudy_sgdGMF_100k_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  script_name <- "CaseStudy_glmpca_100k_avagrad_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  script_name <- "CaseStudy_glmpca_100k_fisher_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  
  script_name <- "CaseStudy_NewWave_100k_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
}


for (i in c(1:5)){
  
  script_name <- "CaseStudy_sgdGMF_200k_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  script_name <- "CaseStudy_glmpca_200k_avagrad_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  script_name <- "CaseStudy_glmpca_200k_fisher_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  
  script_name <- "CaseStudy_NewWave_200k_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
}


for (i in c(1:5)){
  
  script_name <- "CaseStudy_sgdGMF_300k_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  script_name <- "CaseStudy_glmpca_300k_avagrad_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  script_name <- "CaseStudy_glmpca_300k_fisher_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
  
  
  script_name <- "CaseStudy_NewWave_300k_memory.R"
  script_path <- paste(path_casestudy_memory, script_name, sep="/")
  
  
  cmd_logic <- paste0("setwd('", path_casestudy_memory, "'); source('", script_name, "')")
  
  system2("Rscript", 
          args = c("-e", shQuote(cmd_logic), i),
          stdout = "", 
          stderr = "")
}


# These scripts are the old version without memory computation.
#
# source(paste(path_casestudy, "CaseStudy_sgdGMF_100k_1000_1000_250.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_sgdGMF_200k_2000_1000_250.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_sgdGMF_300k_3000_1000_250.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_glmpca_100k_avagrad_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_glmpca_200k_avagrad_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_glmpca_300k_avagrad_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_glmpca_100k_fisher_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_glmpca_200k_fisher_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_glmpca_300k_fisher_def.R", sep="/"))
# 
# 
# source(paste(path_casestudy, "CaseStudy_NewWave_100k_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_100k_def_2.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_100k_def_3.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_100k_def_4.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_100k_def_5.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_200k_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_200k_def_2.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_200k_def_3.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_200k_def_4.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_200k_def_5.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_300k_def.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_300k_def_2.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_300k_def_3.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_300k_def_4.R", sep="/"))
# source(paste(path_casestudy, "CaseStudy_NewWave_300k_def_5.R", sep="/"))


# This script computes the heatmap plot of the manuscript.
source(paste(path_casestudy, "Heatmap_plot.R", sep="/"))

# This script computes the time plot of the manuscript.
source(paste(path_casestudy, "time_plots.R", sep="/"))

# This script computes the memory plot of the manuscript.
source(paste(path_casestudy, "Memory/Memory_plot_peakRAM.R", sep="/"))

# Here, both figures are merged into one figure for the manuscript.
source(paste(path_casestudy, "CaseStudy_merging_plots.R", sep="/"))
