# sgdGMF_Paper
sgdGMF is an R package for efficient estimation of generalized matrix factorization (GMF) models [[1,2,3]](#1,#2,#3).
The package implments the adaptive stochastic gradient descent with block- and coordinate-wise subsampling strategies proposed in [[4]](#4).
Additionally, sgdGMF implements the alternated iterative re-weighted least squares [[1,3]](#1,#3) and diagonal-Hessian quasi-Newton [[1]](#1) algorithms. The package is available at https://github.com/CristianCastiglione/sgdGMF. 

This repository contains a copy of the package and all code required to replicate the results in [[4]](#4).

## Simulation benchmarks

## Arigoni benchmarks

To reproduce the results of the Arigoni dataset, one should do the following steps:

- Download the raw data at: https://figshare.com/articles/dataset/BE1_10XGenomics_count_matrices/23939481/1?file=42312711
- Put all the respective files of each cell-type in Benchmarking/Arigoni/Data/BE1/raw

Run the following R-scripts:
- Benchmarking/Arigoni/preprocessing.R
- Benchmarking/Arigoni/Data_preparation.R
- Benchmarking/Arigoni/Eigenvalues_model_selection.R
- Benchmarking/Arigoni/sgdGMF-cv.R, Benchmarking/Arigoni/sgdGMF-fitting.R
- Benchmarking/Arigoni/NewWave_comparison.R, Benchmarking/Arigoni/Figure_preparation.R
- Benchmarking/Arigoni/Figures_construction.R


## CaseStudy benchmarks
To reproduce the results of the CaseStudy, one should do the following steps:
Run the following R-scripts:
- TENxBrainData_Processing.R
- High_variable_genes.R
- Model_selection_eigenvalues.R
- CaseStudy_Model_fitting_full.R
- CaseStudy_sgdGMF_100k_1000_1000_250.R, CaseStudy_sgdGMF_200k_2000_1000_250.R, CaseStudy_sgdGMF_300k_3000_1000_250.R, 
- CaseStudy_glmpca_100k_avagrad_def.R, CaseStudy_glmpca_200k_avagrad_def.R,
CaseStudy_glmpca_300k_avagrad_def.R
- CaseStudy_glmpca_100k_fisher_def.R, CaseStudy_glmpca_200k_fisher_def.R, CaseStudy_glmpca_300k_fisher_def.R
- CaseStudy_NewWave_100k_def.R, CaseStudy_NewWave_200k_def.R, CaseStudy_NewWave_300k_def.R and the variants followed by _2, _3, _4 and_5
- Heatmap_plot.R
- time_plots.R
- CaseStudy_merging_plots.R

## References
<a id="1">[1]</a>
Collins, M., Dasgupta, S., Schapire, R.E. (2001).
A generalization of principal components analysis to the exponential family.
Advances in neural information processing systems, 14.

<a id="2">[2]</a>
Kidzinski, L., Hui, F.K.C., Warton, D.I., Hastie, T.J. (2022).
Generalized Matrix Factorization: efficient algorithms for fitting generalized linear latent variable models to large data arrays.
Journal of Machine Learning Research, 23(291): 1--29.

<a id="3">[3]</a>
Wang, L., Carvalho, L. (2023).
Deviance matrix factorization.
Electronic Journal of Statistics, 17(2): 3762--3810.

<a id="4">[4]</a>
Castiglione, C., Segers, A., Clement, L, Risso, D. (2024).
Stochastic gradient descent estimation of generalized matrix factorization models with application to single-cell RNA sequencing data.
arXiv preprint: arXiv:2412.20509.

