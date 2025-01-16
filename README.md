# sgdGMF_Paper
sgdGMF is an R package for efficient estimation of generalized matrix factorization (GMF) models [[1,2,3]](#1,#2,#3).
The package implments the adaptive stochastic gradient descent with block- and coordinate-wise subsampling strategies proposed in [[4]](#4).
Additionally, sgdGMF implements the alternated iterative re-weighted least squares [[1,3]](#1,#3) and diagonal-Hessian quasi-Newton [[1]](#1) algorithms. The package is available at https://github.com/CristianCastiglione/sgdGMF. 

This repository contains a copy of the package and all code required to replicate the results in [[4]](#4).

The reproduce the results of the manuscript, one should:
- Download the raw data at: https://figshare.com/articles/dataset/BE1_10XGenomics_count_matrices/23939481/1?file=42312711
- Put all the respective files of each cell-type in Benchmarking/Arigoni/Data/BE1/raw

Then, one can run sgdGMF_manuscript.R script to obtain all results of the paper (this can take a lot of time for some scripts - mainly for the CaseStudy).

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

