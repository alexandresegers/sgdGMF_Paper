
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load our functions
devtools::load_all()
source("sim/utilities.R")

## DATA SIMULATION ----

data.simulation = function (
    n = 1000, m = 100, setting = 1,
    seed = NULL, pca = FALSE, tsne = FALSE
) {

  sim = NULL
  n1 = floor(n/3)
  n2 = floor(n/3)
  n3 = n - n1 - n2
  if (is.null(seed)) { seed = sample.int(100000000, 1) }

  if (setting == 1) { # ELLIPTIC GROUPS
    params = splatter::newSplatParams()
    params = splatter::setParam(params, "batchCells", c(n1, n2, n3))
    params = splatter::setParam(params, "nGenes", m)
    params = splatter::setParam(params, "group.prob", c(0.1, 0.2, 0.2, 0.2, 0.3))
    params = splatter::setParam(params, "de.prob", c(0.3, 0.1, 0.2, 0.01, 0.1))
    params = splatter::setParam(params, "de.downProb", c(0.1, 0.4, 0.9, 0.6, 0.5))
    params = splatter::setParam(params, "de.facLoc", c(0.6, 0.1, 0.1, 0.01, 0.2))
    params = splatter::setParam(params, "de.facScale", c(0.1, 0.4, 0.2, 0.5, 0.4))
    params = splatter::setParam(params, "seed", seed)
    sim = splatter::splatSimulateGroups(params, verbose = FALSE)
    sim = scater::logNormCounts(sim)
  }
  if (setting == 2) { # LINEAR PATHS
    params = splatter::newSplatParams()
    params = splatter::setParam(params, "batchCells", c(n1, n2, n3))
    params = splatter::setParam(params, "nGenes", m)
    params = splatter::setParam(params, "group.prob", c(0.1, 0.2, 0.2, 0.2, 0.3))
    params = splatter::setParam(params, "de.prob", 0.5)
    params = splatter::setParam(params, "de.facLoc", 0.2)
    params = splatter::setParam(params, "path.from", c(0, 1, 2, 3, 4))
    params = splatter::setParam(params, "seed", seed)
    sim = splatter::splatSimulatePaths(params, verbose = FALSE)
    sim = scater::logNormCounts(sim)
  }
  if (setting == 3) { # BRANCHING PATHS
    params = splatter::newSplatParams()
    params = splatter::setParam(params, "batchCells", c(n1, n2, n3))
    params = splatter::setParam(params, "nGenes", m)
    params = splatter::setParam(params, "group.prob", c(0.25, 0.25, 0.25, 0.25))
    params = splatter::setParam(params, "de.prob", 0.5)
    params = splatter::setParam(params, "de.facLoc", 0.2)
    params = splatter::setParam(params, "path.from", c(0, 1, 1, 3))
    params = splatter::setParam(params, "seed", seed)
    sim = splatter::splatSimulatePaths(params, verbose = FALSE)
    sim = scater::logNormCounts(sim)
  }

  logcounts = as.data.frame(logcounts(sim))
  counts = as.data.frame(counts(sim))
  cells = as.data.frame(colData(sim))
  genes = as.data.frame(rowData(sim))
  meta = metadata(sim)

  logpca = NULL
  logtsne = NULL
  if (pca) logpca = RSpectra::svds(scale(t(as.matrix(logcounts))), k = 2)
  if (tsne) logtsne = Rtsne::Rtsne(t(as.matrix(logcounts)), dims = 2, partial_pca = TRUE)

  list(
    logcounts = t(logcounts),
    counts = t(counts),
    cells = cells,
    genes = genes,
    meta = meta,
    logpca = logpca,
    logtsne = logtsne
  )
}

## MAIN DEFINITION ----

main = function (
    niter = 10, setting = 1, nrows = 1000,
    ncols = 100, ncomp = 5, write = FALSE
) {

  # Data dimensions
  n = nrows
  m = ncols
  d = ncomp

  family = poisson()
  error = data.frame()
  verbose = FALSE

  # Define the path and file where to save the results
  # The filename is uniquely identified by the date
  # and time of the main execution
  fileset = NULL
  if (setting == 1) {fileset = "A"}
  if (setting == 2) {fileset = "B"}
  if (setting == 3) {fileset = "C"}
  fileid = format(Sys.time(), "%d-%m-%Y-%H-%M")
  filepath = join.path("sim", "data")
  filename = join.string("summary_", fileset, "_n", n, "_m", m, "_d", d, "_i", niter, ".csv")

  # Simulation loop
  cat(rep("-", 30), sep = "")
  for (iter in 1:niter) {

    # Data simulation
    data = data.simulation(n = n, m = m)
    mask = train.test.split(data$counts, 0.3)
    cells = data$cells$Group

    y = as.matrix(data$counts)
    Z = matrix(1, nrow = m, ncol = 1)
    X = model.matrix(~ Batch, data = data$cells)

    # Train-test split
    train = y
    test = y

    train[mask$train] = NA
    test[mask$test] = NA

    # Matrix completion
    ctrain = naive.completion(train)

    # Model initialization
    cat.alg = function (iter, alg) cat("\n iter:", iter, "  alg:", alg)
    init.alg = function (alg) list (model = alg, u = NULL, v = NULL, d = NULL,
                                    bx = NULL, bz = NULL, eta = NULL, mu = NULL,
                                    tsne = NULL, dev = -1, error = -1, time = NULL)

    # .pearson = init.alg("Pearson")
    # .deviance = init.alg("Deviance")
    .avagrad = init.alg("AvaGrad")
    .fisher = init.alg("Fisher")
    .nbwave = init.alg("NBWaVE")
    .nmf = init.alg("NMF")
    .nnlm = init.alg("NNLM")
    .cmf = init.alg("CMF")
    .airwls = init.alg("AIRWLS")
    .newton = init.alg("Newton")
    .sgd = init.alg("SGD")

    # Model fitting
    cat.alg(iter, "AvaGrad"); try(.avagrad <- fit.avagrad(y=ctrain, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=1000, tol=1e-05))
    cat.alg(iter, "Fisher"); try(.fisher <- fit.fisher(y=ctrain, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=200, tol=1e-05))
    cat.alg(iter, "NBWaVE"); try(.nbwave <- fit.nbwave(y=ctrain, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=500, tol=1e-04))
    cat.alg(iter, "CMF"); try(.cmf <- fit.cmf(y=train, x=X, z=NULL, ncomp=ncomp, family=family, verbose=verbose, maxiter=1000))
    cat.alg(iter, "NMF"); try(.nmf <- fit.nmf(y=ctrain, x=NULL, z=NULL, ncomp=ncomp, family=family, verbose=verbose))
    cat.alg(iter, "NMF+"); try(.nnlm <- fit.nnlm(y=train, x=NULL, z=NULL, ncomp=ncomp, family=family, verbose=verbose, maxiter=2000))
    cat.alg(iter, "AIRWLS"); try(.airwls <- fit.C.airwls(y=train, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=200, stepsize=0.2))
    cat.alg(iter, "Newton"); try(.newton <- fit.C.newton(y=train, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=200, stepsize=0.2))
    cat.alg(iter, "SGD"); try(.sgd <- fit.C.bsgd(y=train, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=500, stepsize=0.01))
    cat("\n", rep("-", 30), sep = "")

    # Model summary
    .err.full  = error.summary(iter, y, cells, .avagrad, .fisher, .nbwave,
                               .cmf, .nmf, .nnlm, .airwls, .newton, .sgd)
    .err.train = error.summary(iter, train, cells, .avagrad, .fisher, .nbwave,
                               .cmf, .nmf, .nnlm, .airwls, .newton, .sgd)
    .err.test  = error.summary(iter, test, cells, .avagrad, .fisher, .nbwave,
                               .cmf, .nmf, .nnlm, .airwls, .newton, .sgd)

    .err.full$dataset  = "Full"
    .err.train$dataset = "Train"
    .err.test$dataset  = "Test"

    .error = rbind(.err.full, .err.train, .err.test)

    # Append the results to the summary error data-frame
    error = rbind(error, .error)

    # Show the partial results
    rownames(error) = 1:nrow(error)
    error$dataset = factor(error$dataset, levels = c("Full", "Train", "Test"))
    error$model = factor(error$model, levels = c(
      "CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "AIRWLS", "Newton", "SGD"))
    #  "AvaGrad", "Fisher", "NBWaVE", "NMF", "NMF+", "CMF", "AIRWLS", "Newton", "SGD"))

    plt = ggplot(data = error, mapping = aes(x = model, color = dataset, fill = dataset)) +
      labs(x = "Models", fill = "", color = "") + theme_bw() + theme(axis.title.x = element_blank())

    plt = ggpubr::ggarrange(
      plt + geom_boxplot(mapping = aes(y = deviance), alpha = 0.5) + labs(y = "Residual deviance"),
      plt + geom_boxplot(mapping = aes(y = error), alpha = 0.5) + labs(y = "Residual squared error"),
      plt + geom_boxplot(mapping = aes(y = silhouette), alpha = 0.5) + labs(y = "Average silhouette"),
      plt + geom_boxplot(mapping = aes(y = log10(time)), alpha = 0.5) + labs(y = "Execution time"),
      nrow = 4, ncol = 1, legend = "right", common.legend = TRUE)

    print(plt)

    # Save all the intermediate results
    if (write) {
      filelist = list.files(path = filepath)
      if (filename %in% filelist) {
        # If the summary file already exist, then just append the summary statistics to the existing file
        write.table(.error, file = join.path(filepath, filename), append = TRUE,
                    col.names = FALSE, row.names = FALSE, sep = ";", dec = ".")
      } else {
        # Otherwise, create a new file and save the summary statistics therein
        write.table(error, file = join.path(filepath, filename), append = FALSE,
                    col.names = TRUE, row.names = FALSE, sep = ";", dec = ".")
      }
    }
  }

  cat("\n")

  return (error)
}

## MAIN EXECUTION ----

RUN = FALSE

if (RUN) main(niter = 25, setting = 1, nrows = 10 *  100, ncol =  100, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 *  250, ncol =  250, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 *  500, ncol =  500, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 *  750, ncol =  750, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 * 1000, ncol = 1000, ncomp = 5, write = TRUE)

if (RUN) main(niter = 25, setting = 2, nrows = 10 *  100, ncol =  100, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 *  250, ncol =  250, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 *  500, ncol =  500, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 *  750, ncol =  750, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 * 1000, ncol = 1000, ncomp = 5, write = TRUE)

if (RUN) main(niter = 25, setting = 3, nrows = 10 *  100, ncol =  100, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 *  250, ncol =  250, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 *  500, ncol =  500, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 *  750, ncol =  750, ncomp = 5, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 * 1000, ncol = 1000, ncomp = 5, write = TRUE)

if (RUN) main(niter = 25, setting = 1, nrows = 10 * 500, ncol = 500, ncomp = 10, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 * 500, ncol = 500, ncomp = 15, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 * 500, ncol = 500, ncomp = 20, write = TRUE)
if (RUN) main(niter = 25, setting = 1, nrows = 10 * 500, ncol = 500, ncomp = 25, write = TRUE)

if (RUN) main(niter = 25, setting = 2, nrows = 10 * 500, ncol = 500, ncomp = 10, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 * 500, ncol = 500, ncomp = 15, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 * 500, ncol = 500, ncomp = 20, write = TRUE)
if (RUN) main(niter = 25, setting = 2, nrows = 10 * 500, ncol = 500, ncomp = 25, write = TRUE)

if (RUN) main(niter = 25, setting = 3, nrows = 10 * 500, ncol = 500, ncomp = 10, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 * 500, ncol = 500, ncomp = 15, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 * 500, ncol = 500, ncomp = 20, write = TRUE)
if (RUN) main(niter = 25, setting = 3, nrows = 10 * 500, ncol = 500, ncomp = 25, write = TRUE)

