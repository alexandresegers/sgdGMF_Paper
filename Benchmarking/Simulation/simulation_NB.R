
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

# Working path
DIRPATH  = "Benchmarking/Simulation"
IMGPATH  = "Benchmarking/Simulation/img"
DATAPATH = "Benchmarking/Simulation/data"

# Global variables
TEST = FALSE
RUN = FALSE

## Load our functions
source(paste(DIRPATH, "utilities.R", sep = "/"))

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

## RUN SIMULATION ----

run.simulation = function (
    niter = 10, from = 0, setting = 1, nrows = 1000,
    ncols = 100, ncomp = 5, write = FALSE
) {

  # Data dimensions
  n = nrows
  m = ncols

  nbfamily = MASS::neg.bin(theta=10)
  family = poisson()
  error = data.frame()
  verbose = FALSE

  # Define the path and file where to save the results
  # The filename is uniquely identified by the date
  # and time of the main execution
  fileset = NULL
  if (setting == 1) {fileset = "bubble"}
  if (setting == 2) {fileset = "linear"}
  if (setting == 3) {fileset = "branch"}
  fileid = format(Sys.time(), "%d-%m-%Y-%H-%M")
  filepath = join.path("Benchmarking", "Simulation", "data", "sim_sgdgmf_nb")
  filename = join.string("summary_sim_n", n, "_m", m, "_d", ncomp, "_i", niter + from, "_NB.csv")

  # Model initialization
  cat.alg = function (iter, alg) cat("\n iter:", iter, "  alg:", alg)
  init.alg = function (alg) list (model = alg, u = NULL, v = NULL, d = NULL,
                                  bx = NULL, bz = NULL, eta = NULL, mu = NULL,
                                  tsne = NULL, dev = -1, error = -1, time = NULL)

  # Simulation loop
  cat(rep("-", 40), "\n", sep = "")
  cat(" setting = ", setting, ", n = ", n, ", m = ", m, ", d = ", ncomp, "\n", sep = "")
  cat(rep("-", 40), sep = "")
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

    # Fitting objects initialization
    .airwls_PO = init.alg("AIRWLS-Pois")
    .newton_PO = init.alg("Newton-Pois")
    .sgd_PO = init.alg("SGD-Pois")
    .airwls_NB = init.alg("AIRWLS-NB")
    .newton_NB = init.alg("Newton-NB")
    .sgd_NB = init.alg("SGD-NB")

    # Model fitting
    cat.alg(iter + from, "AIRWLS-Pois"); try(.airwls_PO <- fit.airwls(y=train, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=200, stepsize=0.1))
    cat.alg(iter + from, "Newton-Pois"); try(.newton_PO <- fit.newton(y=train, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=200, stepsize=0.2))
    cat.alg(iter + from, "SGD-Pois"); try(.sgd_PO <- fit.block.sgd(y=train, x=X, z=Z, ncomp=ncomp, family=family, verbose=verbose, maxiter=500, stepsize=0.01))
    cat.alg(iter + from, "AIRWLS-NB"); try(.airwls_NB <- fit.airwls(y=train, x=X, z=Z, ncomp=ncomp, family=nbfamily, verbose=verbose, maxiter=200, stepsize=0.1))
    cat.alg(iter + from, "Newton-NB"); try(.newton_NB <- fit.newton(y=train, x=X, z=Z, ncomp=ncomp, family=nbfamily, verbose=verbose, maxiter=200, stepsize=0.2))
    cat.alg(iter + from, "SGD-NB"); try(.sgd_NB <- fit.block.sgd(y=train, x=X, z=Z, ncomp=ncomp, family=nbfamily, verbose=verbose, maxiter=500, stepsize=0.01))
    cat("\n", rep("-", 30), sep = "")

    .airwls_PO$model <- "AIRWLS-Pois"
    .newton_PO$model <- "Newton-Pois"
    .sgd_PO$model <- "SGD-Pois"
    .airwls_NB$model <- "AIRWLS-NB"
    .newton_NB$model <- "Newton-NB"
    .sgd_NB$model <- "SGD-NB"

    # Model summary
    .err.full  = error.summary(y, cells, .airwls_PO, .newton_PO, .sgd_PO, .airwls_NB, .newton_NB, .sgd_NB)
    .err.train = error.summary(train, cells, .airwls_PO, .newton_PO, .sgd_PO, .airwls_NB, .newton_NB, .sgd_NB)
    .err.test  = error.summary(test, cells, .airwls_PO, .newton_PO, .sgd_PO, .airwls_NB, .newton_NB, .sgd_NB)

    .err.full$Set  = "Full"
    .err.train$Set = "Train"
    .err.test$Set  = "Test"

    .err.full$Iter  = iter + from
    .err.train$Iter = iter + from
    .err.test$Iter  = iter + from

    .error = rbind(.err.full, .err.train, .err.test)

    # Append the results to the summary error data-frame
    error = rbind(error, .error)

    # Save all the intemediate results
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

    rownames(error) = 1:nrow(error)
    error$Set = factor(error$Set, levels = c("Full", "Train", "Test"))
    error$Model = factor(error$Model, levels = c("AIRWLS-Pois", "Newton-Pois", "SGD-Pois",
                                                 "AIRWLS-NB", "Newton-NB", "SGD-NB"))

    plt = ggplot(data = error, mapping = aes(x = Model, color = Set, fill = Set)) +
      labs(fill = "", color = "") + theme_bw() +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

    plt = ggpubr::ggarrange(
      plt + geom_boxplot(mapping = aes(y = Dev), alpha = 0.5) + labs(title = "Residual deviance"),
      plt + geom_boxplot(mapping = aes(y = RSS), alpha = 0.5) + labs(title = "Residual sum of squares"),
      plt + geom_boxplot(mapping = aes(y = Sil), alpha = 0.5) + labs(title = "Average silhouette"),
      plt + geom_boxplot(mapping = aes(y = Purity), alpha = 0.5) + labs(title = "Average purity"),
      plt + geom_boxplot(mapping = aes(y = log10(Time)), alpha = 0.5) + labs(title = "Execution time"),
      plt + geom_boxplot(mapping = aes(y = log2(Memory)), alpha = 0.5) + labs(title = "Memory consumption"),
      nrow = 3, ncol = 2, legend = "right", common.legend = TRUE, align = "hv")

    print(plt)

    # Free the memory
    rm(data, mask, cells, y, X, Z, test, train, ctrain, .error)
    gc(verbose = FALSE)
  }

  cat("\n")

  return (error)
}

## TEST FUNCTION ----

test.simulation = function () {
  niter = 5
  from = 0
  setting = 1
  write = TRUE
  nrows = 1000
  ncols = 100
  ncomp = 5

  run.simulation(niter = niter, from = from, setting = setting,
                 nrows = nrows, ncols = ncols, ncomp = ncomp, write = write)
}

if (TEST) test.simulation()

## MAIN FUNCTION ----

main = function () {
  niter = 100
  from = 0
  setting = 1
  write = TRUE

  ndata_set_A = c(100, 250, 500, 750, 1000)
  ncomp_set_A = c(5, 5, 5, 5, 5)
  ndata_set_B = c(500, 500, 500, 500)
  ncomp_set_B = c(10, 15, 20, 25)

  # setting A (n and m increasing, d fixed, spherical clusters)
  for (k in 1:5) {
    nrows = ndata_set_A[k] * 10
    ncols = ndata_set_A[k]
    ncomp = ncomp_set_A[k]
    run.simulation(niter = niter, from = from, setting = setting,
                   nrows = nrows, ncols = ncols, ncomp = ncomp, write = write)
  }

  # setting A (n and m fixed, d increasing, spherical clusters)
  for (k in 1:4) {
    nrows = ndata_set_B[k] * 10
    ncols = ndata_set_B[k]
    ncomp = ncomp_set_B[k]
    run.simulation(niter = niter, from = from, setting = setting,
                   nrows = nrows, ncols = ncols, ncomp = ncomp, write = write)
  }
}

if (RUN) main()

## END OF FILE ----

