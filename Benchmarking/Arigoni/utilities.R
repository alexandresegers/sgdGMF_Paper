
## WORKSPACE SETUP ----

# Clean the workspace
# rm(list = ls())
# graphics.off()

# Load the GMF functions
# devtools::load_all()

# Import the needed packages
suppressPackageStartupMessages({
  # Latent variable modelling
  require(MASS)
  #library(gllvm)
  library(glmpca)
  library(NewWave)
  library(COAP)



  # Low-dimensional embedding
  library(Rtsne)
  library(umap)

  # Omics classes and methods
  library(scry)
  library(scater)
  #library(zinbwave)

  # Clustering
  library(cluster)

  # Visualization
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(GGally)
  library(factoextra)
  library(peakRAM)
})

## GLOBAL VARIABLES ----
NCORES = 4

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

rdirichlet = function (n, m = 2, alpha = 1) {
  x = rgamma(n*m, shape = alpha, rate = 1)
  x = matrix(x, nrow = m, ncol = n)
  y = t(apply(x, 2, function(x) x / sum(x)))
  return (y)
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
  # nulldev <- mean(family$dev.resids(y, m, 1), na.rm = TRUE)
  # fitdev <- mean(family$dev.resids(y, x, 1), na.rm = TRUE)
  dev0 <- sgdGMF:::matrix.deviance(m, y, family = family)
  devf <- sgdGMF:::matrix.deviance(x, y, family = family)
  return (devf / dev0)
}

avgsil <- function (y, g) {
  mean(cluster::silhouette(as.numeric(g), dist(y))[,3])
}

# Summary error matrix
error.matrix = function (y, ...) {
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
error.summary = function (y, g, ...) {
  object.list = list(...)
  error = data.frame(Model = c(), Time = c(), Memory = c(), RSS = c(), Cos = c(), Dev = c(), Sil = c())
  for (k in 1:length(object.list)) {
    .object = object.list[[k]]
    .mu = .object$mu
    .model = .object$model
    .time = ifelse(is.null(.object$time), NA, round(.object$time[3], 4))
    .memory = ifelse(is.null(.object$memory), NA, round(unlist(.object$memory[4]),1))
    .rss = ifelse(is.null(.mu), NA, round(rss(y, .mu), 4))
    .cos = ifelse(is.null(.mu), NA, round(cosdist(y, .mu), 4))
    .dev = ifelse(is.null(.mu), NA, round(expdev(y, .mu, family = poisson()), 4))
    .sil = ifelse(is.null(.object$tsne), NA, round(avgsil(.object$tsne, g), 4))
    .error = data.frame(Model = .model, Time = .time, Memory = .memory, RSS = .rss, Cos = .cos, Dev = .dev, Sil = .sil)
    error = rbind(error, .error)
    rownames(error)[k] = k
  }
  return (error)
}

## POSTPROCESSING AND PLOTTING ----

## Plotting function
plot.coeff.matrix <- function (scores, limits = NULL, colours = NULL,
                               transpose = TRUE, symmetric = TRUE) {

  if (is.null(colours)) colours = c("navyblue", "grey95", "darkred")

  if (transpose) scores <- t(scores)

  scores[scores > quantile(scores, .99, na.rm = TRUE)] <- quantile(scores, .99, na.rm = TRUE)
  scores[scores < quantile(scores, .01, na.rm = TRUE)] <- quantile(scores, .01, na.rm = TRUE)

  df <- expand.grid(x = 1:nrow(scores), y = 1:ncol(scores))
  df$z <- as.vector(scores)

  if (is.null(limits)) {
    if (symmetric) {
      limits <- c(-1,+1) * max(abs(scores))
    } else {
      limits <- range(scores)
    }
  }

  plt <- ggplot(data = df, mapping = aes(x = x, y = y, fill = z)) +
    geom_raster() + labs(x = "", y = "", fill = "") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colours, limits = limits) +
    theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank())
  plt
}

# Compare the latent representations
plot.eigvectors = function (u, idx = 1:2, title = "Biplot") {
  ggplot(data = data.frame(x = u[,idx[1]], y = u[,idx[2]], z = 1:nrow(u))) +
    geom_hline(yintercept = 0, lty = 2, color = 1) +
    geom_vline(xintercept = 0, lty = 2, color = 1) +
    geom_point(mapping = aes(x = x, y = y, color = z), size = 2) +
    geom_text(mapping = aes(x = x, y = y - 0.025, color = z, label = z), size = 3) +
    scale_color_gradientn(colors = c("navyblue", "gold"), name = "") +
    labs(x = "PC1", y = "PC2", title = title) + theme_bw() +
    theme(axis.title = element_blank())
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

## Matrix completion for gllvm and glmpca
matrix.completion <- function (y, x = NULL, z = NULL, ncomp = 2,
                               family = poisson(), niter = 10) {

  # data dimensions
  n = nrow(y)
  m = ncol(y)
  d = ncomp
  p = if (is.null(x)) 0 else ncol(x)
  q = if (is.null(z)) 0 else ncol(z)

  # Initialization via iterated least squares
  init = sgdGMF:::ls.svd.init(Y = y, X = x, Z = z, d = d, family = family)

  # create proxy model matrices if they are null
  if (p == 0) x = matrix(0, nrow = n, ncol = p)
  if (q == 0) z = matrix(0, nrow = q, ncol = m)

  # compute the linear predictor
  eta = matrix(0, nrow = n, ncol = m)
  eta[] = eta + tcrossprod(init$u, init$v)
  eta[] = eta + tcrossprod(x, init$bx)
  eta[] = eta + tcrossprod(init$bz, z)

  # compute the estimated mean, variance and deviance residuals
  mu = family$linkinv(eta)
  var = family$variance(mu)

  # fill the missing values
  isna = is.na(y)
  y[isna] = round(mu[isna])

  # compute the deviance and Pearson residuals
  dr = family$dev.resids(y, mu, 1)
  pr = (y - mu) / sqrt(var)

  # return the output
  # list(n = n, m = m, p = p, q = q,
  #      u = init$u, v = init$v,
  #      beta.x = init$bx, beta.z = init$bz,
  #      y = y, eta = eta, mu = mu, var = var,
  #      dr = dr, pr = pr, isna = isna)

  # return the output
  return (y)
}

naive.completion = function (y) {
  apply(y, 2, function (x) {
    m = mean(x, na.rm = TRUE)
    x[is.na(x)] = round(m)
    return (x)
  })
}

## PREPROCESSING ----

null.residuals = function (
    y, x = NULL, z = NULL,
    family = poisson(), type = "deviance"
) {
  n = nrow(y); m = ncol(y)
  y = as.matrix(y)
  yr = rowMeans(y) # row-average
  yc = colMeans(y) # col-average
  yt = mean(yr) # total average
  mu = tcrossprod(yr, yc) / yt

  res = matrix(NA, nrow = n, ncol = m)
  if (type == "deviance") {
    res[] = sign(y - mu) * sqrt(family$dev.resids(y, mu, 1))
  }
  if (type == "pearson") {
    res[] = (y - mu) / sqrt(family$variance(mu))
  }
  if (!is.null(x)) {
    beta = matrix(NA, nrow = m, ncol = ncol(x))
    beta[] = t(solve(crossprod(x, x), crossprod(x, res)))
    res[] = res - tcrossprod(x, beta)
  }
  if (!is.null(z)) {
    gamma = matrix(NA, nrow = n, ncol = ncol(z))
    gamma[] = t(solve(crossprod(z, z), crossprod(z, t(res))))
    res[] = res - tcrossprod(gamma, z)
  }
  return (res)
}

## MODEL FIT ----



# GLMPCA via AVAGRAD
fit.glmpca = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    method = c("avagrad", "fisher"),
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
    optimizer = method,
    ctl = list(verbose = verbose,
               maxIter = maxiter,
               tol = tol)) %>%
    suppressWarnings())
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

  # Explained deviance metrics
  # rss  =  RSS(y, mu)
  # rmse = RMSE(y, mu)
  # cosd = COSD(y, mu)
  # dev  = RDEV(y, mu, family = family)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = switch(method, "avagrad" = "AvaGrad", "fisher" = "Fisher"),
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coefZ,
    bz = fit$coefX,
    eta = eta,
    mu = mu,
    tsne = tsne$Y,
    dev = fit$dev,
    error = error,
    time = timef - time0,
    memory = memory)
}

## NBWaVE
fit.nbwave = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    maxiter = 100, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL,
    NCORES = 4) {

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
    suppressWarnings())
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

  # Explained deviance metrics
  # rss  =  RSS(y, mu)
  # rmse = RMSE(y, mu)
  # cosd = COSD(y, mu)
  # dev  = RDEV(y, mu, family = family)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

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
    tsne = tsne$Y,
    dev = NULL,
    error = error,
    time = timef - time0,
    memory = memory)
}


fit.coap = function(
    y, x = NULL, z = NULL, ncomp = 2,
    maxiter = 1000, tol = 1e-5,
    verbose = FALSE, ncores = 4) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

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
      fast_svd = FALSE) %>%
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
  mu = poisson()$linkinv(eta)

  # Orthogonalization
  uv = svd::propack.svd(uv, neig = ncomp)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = ncores, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "COAP",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    beta = b,
    alpha = a,
    eta = eta,
    mu = mu,
    tsne = scale(tsne$Y),
    dev = NULL,
    error = error,
    time = timef - time0,
    memory = memory)
}

## Averaged stochastic gradient descent
fit.R.bsgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 5000, stepsize = 0.001, tol = 1e-04,
    verbose = FALSE, train = NULL, test = NULL) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  # model fitting
  time0 = proc.time()
  fit = sgdgmf(
    Y = y, X = x, Z = z,
    ncomp = ncomp,
    family = family,
    method = "b-sgd",
    penalty = list(u = penalty, v = penalty, b = 0),
    init = list(init = "svd", niter = 10),
    control = list(maxiter = maxiter,
                   size = c(100, 20),
                   rate0 = stepsize,
                   verbose = verbose,
                   frequency = 250)) %>%
    suppressWarnings()
  timef = proc.time()

  # Latent factor normalization
  uv = tcrossprod(fit$coef$U, fit$coef$V)
  uv = svd::propack.svd(uv, neig = ncomp)

  # Fitted values
  eta = fit$pred$eta
  mu = fit$pred$mu

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c( RSS(train, mu),  RSS(test, mu)),
      cosd = c(COSD(train, mu), COSD(test, mu)),
      dev  = c(RDEV(train, mu, family = family),
               RDEV(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "B-SGD",
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = fit$coef$betaX,
    bz = fit$coef$betaZ,
    eta = fit$pred$eta,
    mu = fit$pred$mu,
    dev = fit$trace$deviance,
    error = error,
    time = timef - time0)
}


## MODEL FIT (C++) ----


## GMF via block-wise SGD
fit.C.bsgd = function (
    y, x = NULL, z = NULL, ncomp = 2, family = poisson(),
    penalty = 1, maxiter = 200, stepsize = 0.01, tol = 1e-05,
    verbose = FALSE, train = NULL, test = NULL, size1 = 100, size2 = 100,
    savedata = FALSE) {

  n = nrow(y)
  m = ncol(y)

  if (is.null(x)) x = matrix(1, nrow = n, ncol = 1)
  if (is.null(z)) z = matrix(1, nrow = m, ncol = 1)

  p = ncol(x)
  q = ncol(z)

  # init = gmf.init(y, x, z, d = ncomp, method = "svd",
  #                 niter = 0, verbose = FALSE)

  familyname = family$family
  linkname = family$link

  if (familyname == "Negative Binomial") {
    familyname = "negbinom"
  }
  time0 = proc.time()
  memory = peakRAM::peakRAM(
    fit <- sgdGMF::sgdgmf.fit(Y = y, X = x, Z = z, ncomp = ncomp, family = family,
                           control.alg = list(maxiter = maxiter, size = c(size1, size2),
                                              savedata = savedata, normalize = F),
                           control.init = list(method = "light", type = "link",
                                               normalize = F),
                           method = "sgd", sampling = "block")
    )
  timef = proc.time()

  # init = sgdGMF::init.param(
  #   Y = y, X = x, Z = z, ncomp = ncomp, family = family,
  #   method = "ols", type = "link", verbose = FALSE)
  #
  # # model fitting
  # fit = sgdGMF::cpp.fit.bsgd(
  #   Y = y, X = x, B = init$B, A = init$A, Z = z, U = init$U, V = init$V,
  #   ncomp = ncomp, familyname = familyname, linkname = linkname,
  #   lambda = c(0,0,1,0), maxiter = maxiter, rate0 = stepsize,
  #   size1 = size1, size2 = size2, eps = 1e-08, nafill = 1, tol = tol,
  #   damping = 1e-03, verbose = verbose, frequency = 100)

  # bz = fit$U[, (p+1):(p+q)]
  # bx = fit$V[, 1:p]
  # u = fit$U[, (p+q+1):(p+q+ncomp)]
  # v = fit$V[, (p+q+1):(p+q+ncomp)]

  bz = fit$A
  bx = fit$B
  u = fit$U
  v = fit$V

  # Latent factor normalization
  uv = tcrossprod(u, v)
  uv = svd::propack.svd(uv, neig = ncomp)

  eta <- x %*% t(fit$B) + fit$A %*% t(z) + fit$U %*% t(fit$V)
  mu <- exp(eta)

  # tSNE low dimensional embedding
  tsne = Rtsne::Rtsne(uv$u, dims = 2, partial_pca = FALSE,
                      num_threads = NCORES, verbose = verbose)

  error = NULL
  if (!is.null(train) && !is.null(test)) {
    error = data.frame(
      rss  = c(rss(train, mu), rss(test, mu)),
      cosd = c(cosdist(train, mu), cosdist(test, mu)),
      dev  = c(expdev(train, mu, family = family),
               expdev(test, mu, family = family)))

    colnames(error) = c("RSS", "Cos", "Dev")
    rownames(error) = c("Train", "Test")
  }

  # Output
  list(
    model = "B-SGD",
    fit = fit,
    u = uv$u,
    v = uv$v,
    d = uv$d,
    bx = bx,
    bz = bz,
    eta = eta,
    mu = mu,
    phi = fit$phi,
    tsne = tsne$Y,
    dev = fit$trace[,2],
    error = error,
    time = timef - time0,
    memory = memory)
}




