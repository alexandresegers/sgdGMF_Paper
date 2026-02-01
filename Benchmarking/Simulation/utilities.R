
## WORKSPACE SETUP ----

# Clean the workspace
# rm(list = ls())
# graphics.off()

# Import the needed packages
suppressPackageStartupMessages({
  # Latent variable modelling
  require(MASS)

  # Fast linear algebra
  library(irlba)
  library(svd)
  library(RSpectra)

  # sgdGMF
  library(sgdGMF)

  # Genalized PCA
  library(gmf)
  library(glmpca)
  library(PoissonPCA)

  # Generalized factor models
  library(gllvm)
  library(GFM)
  library(COAP)
  library(rrpack)

  # Matrix factorization
  library(NMF)
  library(NNLM)
  library(cmfrec)

  # Low-dimensional embedding
  library(Rtsne)
  library(umap)

  # Omics classes and methods
  library(scry)
  library(splatter)
  library(scater)
  library(zinbwave)
  library(NewWave)

  # Clustering
  library(cluster)
  library(bluster)

  # Visualization
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(GGally)
  library(factoextra)
})

## GLOBAL VARIABLES ----
NCORES = 6

## UTILITIES ----

## Join several paths
join.path <- function (...) {
  paste(..., sep = "/")
}

## Join several strings
join.string <- function (...) {
  paste(..., sep = "")
}

## TRANSFORMATIONS ----

minmax <- function (x, a = -1, b = +1) {
  minx <- min(x); maxx <- max(x)
  (b - a) * (x - minx) / (maxx - minx) + a
}

boxcox <- function (x, lambda) {
  (x^lambda - 1) / lambda
}

log1p <- function (x) {
  log(1 + x)
}

log10p1 <- function (x) {
  log10(1 + x)
}

slog10p1 <- function (x) {
  sign(x) * log10(1 + abs(x))
}

ssqrt <- function (x) {
  sign(x) * sqrt(abs(x) + 1e-03)
}

## ERROR MEASURES ----

# Residual sum of squares
rss <- function (y, x, f = log1p) {
  m = mean(y, na.rm = TRUE)
  err0 = mean(c(f(y) - f(m))^2, na.rm = TRUE)
  errf = mean(c(f(y) - f(x))^2, na.rm = TRUE)
  return (errf / err0)
}

# Root mean squared error
rmse <- function (y, x, f = log1p) {
  sqrt(mean((f(y) - f(x))^2, na.rm = TRUE))
}

# Cosine distance
cosdist <- function (y, x, f = log1p) {
  xy <- sum(f(y) * f(x), na.rm = TRUE)
  xx <- sqrt(sum(f(y)^2, na.rm = TRUE))
  yy <- sqrt(sum(f(x)^2, na.rm = TRUE))
  return(1 - abs(xy) / (xx * yy))
}

# Explained deviance
expdev <- function (y, x, family) {
  m <- mean(y, na.rm = TRUE)
  dev0 <- sgdGMF:::matrix.deviance(m, y, family = family)
  devf <- sgdGMF:::matrix.deviance(x, y, family = family)
  return (devf / dev0)
}

# Average silhouette
avgsil <- function (y, g) {
  mean(cluster::silhouette(as.numeric(g), dist(y))[,3])
}

# Neighbor purity
purity <- function (u, g) {
  mean(bluster::neighborPurity(u, clusters = g)$purity)
}

# Summary error matrix
error.matrix <- function (y, ...) {
  object.list = list(...)
  error = data.frame(Model = c(), Time = c(), RSS = c(), Cos = c(), Dev = c())
  for (k in 1:length(object.list)) {
    .object = object.list[[k]]
    .mu = .object$mu
    .model = .object$model
    .time = round(.object$time[3], 4)
    .rss = round(rss(y, .mu), 4)
    .cos = round(cosdist(y, .mu), 4)
    .dev = round(expdev(y, .mu, family = poisson()), 4)
    .error = data.frame(Model = .model, Time = .time, RSS = .rss, Cos = .cos, Dev = .dev)
    error = rbind(error, .error)
    rownames(error)[k] = k
  }
  return (error)
}

# Summary error matrix + silhouette
error.summary <- function (y, g, ...) {
  object.list = list(...)
  error = data.frame(Model = c(), Time = c(), RSS = c(), Cos = c(), Dev = c(), Sil = c())
  for (k in 1:length(object.list)) {
    .object = object.list[[k]]
    .mu = .object$mu
    .model = .object$model
    .time = ifelse(is.null(.object$time), NA, round(.object$time[3], 4))
    .memory = ifelse(is.null(.object$memory), NA, round(.object$memory$Peak_RAM_Used_MiB[1]))
    .rss = ifelse(is.null(.mu), NA, round(rss(y, .mu), 4))
    .cos = ifelse(is.null(.mu), NA, round(cosdist(y, .mu), 4))
    .dev = ifelse(is.null(.mu), NA, round(expdev(y, .mu, family = poisson()), 4))
    .sil = ifelse(is.null(.object$tsne), NA, round(avgsil(.object$tsne, g), 4))
    .purity = ifelse(is.null(.object$u), NA, round(purity(.object$u, g), 4))
    .error = data.frame(Model = .model, Time = .time, Memory = .memory, RSS = .rss,
                        Cos = .cos, Dev = .dev, Sil = .sil, Purity = .purity)
    error = rbind(error, .error)
    rownames(error)[k] = k
  }
  return (error)
}

## TRAIN-TEST SPLIT ----

## Train-test split
train.test.split <- function (y, test = 0.3) {

  n = nrow(y)
  m = ncol(y)

  idx = cbind(x = sample.int(n = n, size = floor(test * n * m), replace = TRUE),
              y = sample.int(n = m, size = floor(test * n * m), replace = TRUE))

  mask = matrix(0, nrow = n, ncol = m)
  mask[idx] = 1

  train = (mask == 1)
  test = (mask != 1)

  list(train = train, test = test)
}

## Naive matrix completion with column mean imputation
naive.completion <- function (y) {
  apply(y, 2, function (x) {
    m = mean(x, na.rm = TRUE)
    x[is.na(x)] = round(m)
    return (x)
  })
}

## MODEL FIT ----

# Latent feature extraction via Pearson residuals
fit.pearson = function (y, x = NULL, z = NULL, ncomp = 2,
                        family = poisson(), verbose = FALSE) {

  # model fitting
  time0 = proc.time()
  res = scry::nullResiduals(object = as.matrix(y), fam = "poisson", type = "pearson")
  SVD = svd::propack.svd(res, neig = ncomp)
  timef = proc.time()

  if (!is.null(X)) beta.x = matrix(0, nrow = ncol(y), ncol = ncol(x))
  if (!is.null(z)) beta.z = matrix(0, nrow = nrow(y), ncol = ncol(z))

  eta = tcrossprod(SVD$u, SVD$v %*% diag(SVD$d))
  mu = family$linkinv(eta)
  tsne = Rtsne::Rtsne(SVD$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "Pearson",
    u = SVD$u %*% diag(sqrt(SVD$d)),
    v = SVD$v %*% diag(sqrt(SVD$d)),
    d = SVD$d,
    bx = beta.x,
    bz = beta.z,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = -1,
    error = -1,
    time = timef - time0)
}

# Latent feature extraction via deviance residuals
fit.deviance = function (y, x = NULL, z = NULL, ncomp = 2,
                         family = poisson(), verbose = FALSE) {

  # model fitting
  time0 = proc.time()
  res = scry::nullResiduals(object = as.matrix(y), fam = "poisson", type = "deviance")
  SVD = svd::propack.svd(res, neig = ncomp)
  timef = proc.time()

  if (!is.null(X)) beta.x = matrix(0, nrow = ncol(y), ncol = ncol(x))
  if (!is.null(z)) beta.z = matrix(0, nrow = nrow(y), ncol = ncol(z))

  eta = tcrossprod(SVD$u, SVD$v %*% diag(SVD$d))
  mu = family$linkinv(eta)
  tsne = Rtsne::Rtsne(SVD$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "Deviance",
    u = SVD$u %*% diag(sqrt(SVD$d)),
    v = SVD$v %*% diag(sqrt(SVD$d)),
    d = SVD$d,
    bx = beta.x,
    bz = beta.z,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = -1,
    error = -1,
    time = timef - time0)
}

# GLLVM via variational approximations
fit.gllvm = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-03,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- gllvm::gllvm(
      y = y, X = x, Z = z,
      formula = ~ .,
      num.lv = ncomp,
      family = family,
      method = "EVA",
      reltol = tol) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$lvs, fit$params$theta)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = residuals(fit)$linpred
  mu = family$linkinv(eta)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)


  # Output
  list(
    model = "gllvm",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    error = NULL,
    time = timef - time0,
    memory = memory)
}

# GFM via alternated maximization
fit.gfmam = function(
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  idx = which(apply(y, 2, sd) == 0)
  if (length(idx) > 0) {
    for (j in idx) {
      y[sample(1:n, 1),j] = 1
    }
  }

  # GFM-AM model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- GFM::gfm(
      XList = list(y),
      types = "poisson",
      q = ncomp,
      offset = FALSE,
      dc_eps = tol,
      maxIter = maxiter,
      verbose = verbose,
      algorithm = "AM") %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  b = fit$hmu
  xb = tcrossprod(rep(1, n), b)

  a = rep(0, m)
  az = matrix(0, nrow = n, ncol = m)

  uv = tcrossprod(fit$hH, fit$hB)

  # Fitted values
  eta = xb + az + uv
  mu = family$linkinv(eta)

  # Orthogonalization
  uv = svd::propack.svd(uv, neig = ncomp)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "GFM-AM",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = b,
    bz = a,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

# GFM via variational EM
fit.gfmvem = function(
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  idx = which(apply(y, 2, sd) == 0)
  if (length(idx) > 0) {
    for (j in idx) {
      y[sample(1:n, 1),j] = 1
    }
  }

  # GFM-VEM model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- GFM::gfm(
      XList = list(y),
      types = "poisson",
      q = ncomp,
      offset = FALSE,
      dc_eps = tol,
      maxIter = maxiter,
      verbose = verbose,
      algorithm = "VEM") %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  b = fit$hmu
  xb = tcrossprod(rep(1, n), b)

  a = rep(0, m)
  az = matrix(0, nrow = n, ncol = m)

  uv = tcrossprod(fit$hH, fit$hB)

  # Fitted values
  eta = xb + az + uv
  mu = family$linkinv(eta)

  # Orthogonalization
  uv = svd::propack.svd(uv, neig = ncomp)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "GFM-VEM",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = b,
    bz = a,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

# COAP (covariate augmented Poisson factor model)
fit.coapf = function(
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  idx = which(apply(y, 2, sd) == 0)
  if (length(idx) > 0) {
    for (j in idx) {
      y[sample(1:n, 1),j] = 1
    }
  }

  # COAP model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- COAP::RR_COAP(
      X_count = y,
      Z = x,
      q = ncomp,
      epsELBO = tol,
      maxIter = maxiter,
      verbose = verbose,
      joint_opt_beta = FALSE,
      fast_svd = TRUE) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  b = fit$bbeta
  xb = tcrossprod(x, b)

  a = rep(0, m)
  az = matrix(0, nrow = n, ncol = m)

  uv = tcrossprod(fit$H, fit$B)

  # Fitted values
  eta = xb + az + uv
  mu = family$linkinv(eta)

  # Orthogonalization
  uv = svd::propack.svd(uv, neig = ncomp)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "COAP",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = b,
    bz = a,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

# GLMPCA via AVAGRAD
fit.glmpca = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-06,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # glmPCA model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- glmpca::glmpca(
      Y = y, X = z, Z = x,
      L = ncomp,
      fam = "poi",
      minibatch = "none",
      optimizer = "fisher",
      ctl = list(verbose = verbose,
                 maxIter = maxiter,
                 tol = tol)) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(as.matrix(fit$loadings), as.matrix(fit$factors))
  uv = svd::propack.svd(uv, neig = ncomp)

  uv = list(d = rep(1, ncomp),
            u = as.matrix(fit$loadings),
            v = as.matrix(fit$factors))

  # Fitted values
  mu = predict(fit)
  eta = family$linkfun(mu)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "glmPCA",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = fit$dev,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

# GLMPCA via AVAGRAD
fit.avagrad = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 1000, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # glmPCA model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- glmpca::glmpca(
      Y = y, X = z, Z = x,
      L = ncomp,
      fam = "poi",
      minibatch = "none",
      optimizer = "avagrad",
      ctl = list(verbose = verbose,
                 maxIter = maxiter,
                 tol = tol)) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(as.matrix(fit$loadings), as.matrix(fit$factors))
  uv = svd::propack.svd(uv, neig = ncomp)

  uv = list(d = rep(1, ncomp),
            u = as.matrix(fit$loadings),
            v = as.matrix(fit$factors))

  # Fitted values
  mu = predict(fit)
  eta = family$linkfun(mu)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "AvaGrad",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = fit$dev,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

# GLMPCA via Fisher scoring
fit.fisher = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 100, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # glmPCA model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- glmpca::glmpca(
      Y = y, X = z, Z = x,
      L = ncomp,
      fam = "poi",
      minibatch = "none",
      optimizer = "fisher",
      ctl = list(verbose = verbose,
                 maxIter = maxiter,
                 tol = tol)) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(as.matrix(fit$loadings), as.matrix(fit$factors))
  uv = svd::propack.svd(uv, neig = ncomp)

  uv = list(d = rep(1, ncomp),
            u = as.matrix(fit$loadings),
            v = as.matrix(fit$factors))

  # Fitted values
  mu = predict(fit)
  eta = family$linkfun(mu)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "Fisher",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = fit$dev,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

## NBWaVE
fit.nbwave = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 100, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  idx = which(apply(y, 2, sd) == 0)
  if (length(idx) > 0) {
    for (j in idx) {
      y[sample(1:n, 1),j] = 1
    }
  }

  # NBWaVE model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- NewWave::newFit(
      Y = t(y),
      X = x,
      V = z,
      K = ncomp,
      commondispersion = TRUE,
      verbose = verbose,
      maxiter_optimize = maxiter,
      stop_epsilon = tol,
      children = NCORES) %>%
      suppressWarnings() %>%
      suppressMessages())
  timef = proc.time()

  # Latent factor normalization
  b = t(fit@beta)
  xb = tcrossprod(x, b)

  a = t(fit@gamma)
  az = tcrossprod(a, z)

  uv = fit@W %*% fit@alpha

  # Fitted values
  eta = xb + az + uv
  mu = family$linkinv(eta)

  # Orthogonalization
  uv = svd::propack.svd(uv, neig = ncomp)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "NBWaVE",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = b,
    bz = a,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

## NMF
fit.nmf = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    verbose = FALSE, train = NULL, test = NULL) {

  suppressPackageStartupMessages(require("NMF"))

  time0 = proc.time()
  if (verbose) {
    memory = peakRAM::peakRAM(
      fit <- NMF::nmf(x = y, rank = ncomp, method = "brunet", seed = "nndsvd", nrun = 1, .options = "v") %>%
        suppressWarnings() %>% suppressMessages()
    )
  } else {
    memory = peakRAM::peakRAM(
      fit <- NMF::nmf(x = y, rank = ncomp, method = "brunet", seed = "nndsvd", nrun = 1) %>%
        suppressWarnings() %>% suppressMessages()
    )

  }
  timef = proc.time()

  mu = fit@fit@W %*% fit@fit@H

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(fit@fit@W, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "NMF",
    u = fit@fit@W,
    v = t(fit@fit@H),
    d = NULL,
    bx = NULL,
    bz = NULL,
    eta = NULL,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

## NNLM
fit.nnlm = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 1000, verbose = FALSE,
    train = NULL, test = NULL) {

  suppressPackageStartupMessages(require("NNLM"))

  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- NNLM::nnmf(
      A = y, k = ncomp, alpha = penalty, beta = penalty, n.threads = NCORES,
      method = "lee", loss = "mkl", max.iter = maxiter, verbose = verbose)) %>%
    suppressWarnings() %>% suppressMessages()
  timef = proc.time()

  mu = fit$W %*% fit$H

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(fit$W, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "NMF+",
    u = fit$W,
    v = t(fit$H),
    d = NULL,
    bx = NULL,
    bz = NULL,
    eta = NULL,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = fit$nkl,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

## NNLM
fit.cmf = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 1000, verbose = FALSE,
    train = NULL, test = NULL) {

  suppressPackageStartupMessages(require("cmfrec"))

  if (!is.null(x)) if (sd(x[,1]) == 0) x = x[, -1, drop = FALSE]
  if (!is.null(z)) if (sd(z[,1]) == 0) z = z[, -1, drop = FALSE]

  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- cmfrec::CMF(
      X = y, U = x, I = z, k = ncomp, nonneg = TRUE,
      user_bias = FALSE, item_bias = FALSE, center = FALSE,
      nthreads = NCORES, niter = maxiter, verbose = verbose) %>%
      suppressWarnings() %>% suppressMessages())
  timef = proc.time()

  eta = crossprod(fit$matrices$A, fit$matrices$B) + fit$matrices$glob_mean
  mu = eta

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(t(fit$matrices$A), dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "CMF",
    u = t(fit$matrices$A),
    v = fit$matrices$B,
    d = NULL,
    bx = NULL,
    bz = NULL,
    eta = NULL,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = NULL,
    time = timef - time0,
    memory = memory)
}

## MODEL sgdGMF FIT ----

## GMF via AIRWLS
fit.airwls = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

  # model set-up
  lambda = list(B = 1e-4, A = 1e-4, U = penalty, V = 1e-4)
  ctr.init = list(method = "ols", type = "link")
  ctr.alg = list(normalize = FALSE, maxiter = maxiter, nstep = 1, stepsize = stepsize,
                 eps = 1e-16, nafill = 1, tol = tol,damping = 1e-3, verbos = verbose,
                 frequency = 25, parallel = TRUE, nthreads = NCORES)

  # model fitting
  time0 = proc.time()
  lambda = (1 - 1e-04) * c(0, 0, penalty, 0) + 1e-04
  memory = peakRAM::peakRAM(
    fit <- sgdGMF::sgdgmf.fit(Y = y, X = x, Z = z, family = family,
                              ncomp = ncomp, method = "airwls", penalty = lambda,
                              control.init = ctr.init, control.alg = ctr.alg))
  timef = proc.time()

  bz = fit$A
  bx = fit$B
  u = fit$U
  v = fit$V

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "AIRWLS",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = scale(tsne$Y),
    dev = fit$trace[,2],
    error = NULL,
    time = timef - time0,
    memory = memory)
}


## GMF via quasi-Newton
fit.newton = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

  # model set-up
  lambda = list(B = 0, A = 0, U = penalty, V = 0)
  ctr.init = list(method = "ols", type = "link")
  ctr.alg = list(normalize = FALSE, maxiter = maxiter, stepsize = stepsize,
                 eps = 1e-16, nafill = 1, tol = tol, damping = 1e-3, verbos = verbose,
                 frequency = 25, parallel = TRUE, nthreads = NCORES)

  # model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- sgdGMF::sgdgmf.fit(Y = y, X = x, Z = z, family = family,
                              ncomp = ncomp, method = "newton", penalty = lambda,
                              control.init = ctr.init, control.alg = ctr.alg))
  timef = proc.time()

  bz = fit$A
  bx = fit$B
  u = fit$U
  v = fit$V

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "Newton",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = scale(tsne$Y),
    dev = fit$trace[,2],
    error = NULL,
    time = timef - time0,
    memory = memory)
}


## GMF via block-wise SGD
fit.block.sgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }

  # model set-up
  lambda = list(B = 0, A = 0, U = penalty, V = 0)
  ctr.init = list(method = "ols", type = "link")
  ctr.alg = list(normalize = FALSE, maxiter = maxiter, size = c(100, 20),
                 rate0 = stepsize, eps = 1e-16, nafill = 1, tol = tol,
                 damping = 1e-3, verbos = verbose, frequency = 100)

  # model fitting
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- sgdGMF::sgdgmf.fit(Y = y, X = x, Z = z, family = family, ncomp = ncomp,
                              method = "sgd", sampling = "block", penalty = lambda,
                              control.init = ctr.init, control.alg = ctr.alg))
  timef = proc.time()

  bz = fit$A
  bx = fit$B
  u = fit$U
  v = fit$V

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  mu = fit$mu

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  # Output
  list(
    model = "SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = fit$eta,
    mu = fit$mu,
    phi = fit$phi,
    tsne = scale(tsne$Y),
    dev = fit$trace[,2],
    error = NULL,
    time = timef - time0,
    memory = memory)
}


## END OF FILE ----
