# ----------- Simulation benchmark

# First, all data is simulated with splatter, and all methods are run on these
# simulated datasets.
source("Benchmarking/Simulation/simulation_MAIN.R")
source("Benchmarking/Simulation/simulation_NB.R")

# Then, figures are made for the summary statistics and the tSNE projections.
source("Benchmarking/Simulation/plot_sim_example_MAIN.R")

# Similarly, figures are made for the summary statistics in function of the
# matrix dimensions and latent space rank.
source("Benchmarking/Simulation/plot_sim_summary_MAIN.R")
source("Benchmarking/Simulation/plot_sim_summary_NB.R")


# ----------- Arigoni benchmark

# First, the raw data of all celltypes is merged together in one file, and
# an initial quality control is done.
source("Benchmarking/Arigoni/preprocessing.R")

# Then, the data is additionally filtered. For example, PBMC cells are not used
# and cells with many mitochondrial reads are removed.
source("Benchmarking/Arigoni/Data_preparation.R")

# We calculate a scree-plot based on the eigenvalues for the model selection.
source("Benchmarking/Arigoni/Eigenvalues_model_selection.R")

# We also perform cross-validation based on model selection criteria and
# out of sample deviances to select the number of latent factors.
source("Benchmarking/Arigoni/sgdGMF-cv.R")

# sgdGMF estimation is performed using the optimal number of latent factors
# chosen in the cross-validation.
source("Benchmarking/Arigoni/sgdGMF-fitting.R")

# Benchmarking with glmPCA and NewWave is done. Here, these models are
# used to perform dimensionality reduction on the Arigoni dataset.
source("Benchmarking/Arigoni/NewWave_comparison.R")

# We here group all results in data.frames that can be used to make the plots
# of the manuscript.
source("Benchmarking/Arigoni/Figure_preparation.R")

# We here construct the plots of the manuscript.
source("Benchmarking/Arigoni/Figures_construction.R")



# ----------- CaseStudy benchmark

# First, the dataset is downloaded and quality control is performed to
# remove aberrant samples or genes.
source("Benchmarking/CaseStudy/TENxBrainData_Processing.R")

# The variance of each gene is estimated to select the most variable
# genes of the dataset. These genes will be further used in the decomposition.
source("Benchmarking/CaseStudy/High_variable_genes.R")

# We calculate a scree-plot based on the eigenvalues for the model selection.
source("Benchmarking/CaseStudy/Model_selection_eigenvalues.R")

# sgdGMF estimation is performed using the optimal number of latent factors
# estimated with the scree plot.
source("Benchmarking/CaseStudy/CaseStudy_Model_fitting_full.R")

# For the comparisons of the time-usage, sgdGMF is also run on smaller
# parts of the data. Here, 5 times, the time is measured to compute the
# dimensionality reduction for 100.000, 200.000 and 300.000 cells in this
# dataset. This is also done for glmPCA with Avagrad estimation and with Fisher
# estimation. Also, NewWave is used on the subsampled dataset. However, as
# NewWave couldn't be run 5 times within the limited computing time (72 hours),
# it was run in 5 different scripts for each subsampled dataset.

source("Benchmarking/CaseStudy/CaseStudy_sgdGMF_100k_1000_1000_250.R")
source("Benchmarking/CaseStudy/CaseStudy_sgdGMF_200k_2000_1000_250.R")
source("Benchmarking/CaseStudy/CaseStudy_sgdGMF_300k_3000_1000_250.R")
source("Benchmarking/CaseStudy/CaseStudy_glmpca_100k_avagrad_def.R")
source("Benchmarking/CaseStudy/CaseStudy_glmpca_200k_avagrad_def.R")
source("Benchmarking/CaseStudy/CaseStudy_glmpca_300k_avagrad_def.R")
source("Benchmarking/CaseStudy/CaseStudy_glmpca_100k_fisher_def.R")
source("Benchmarking/CaseStudy/CaseStudy_glmpca_200k_fisher_def.R")
source("Benchmarking/CaseStudy/CaseStudy_glmpca_300k_fisher_def.R")


source("Benchmarking/CaseStudy/CaseStudy_NewWave_100k_def.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_100k_def_2.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_100k_def_3.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_100k_def_4.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_100k_def_5.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_200k_def.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_200k_def_2.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_200k_def_3.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_200k_def_4.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_200k_def_5.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_300k_def.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_300k_def_2.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_300k_def_3.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_300k_def_4.R")
source("Benchmarking/CaseStudy/CaseStudy_NewWave_300k_def_5.R")


# This script computes the heatmap plot of the manuscript.
source("Benchmarking/CaseStudy/Heatmap_plot.R")

# This script computes the time plot of the manuscript.
source("Benchmarking/CaseStudy/time_plots.R")

# Here, both figures are merged into one figure for the manuscript.
source("Benchmarking/CaseStudy/CaseStudy_merging_plots.R")
