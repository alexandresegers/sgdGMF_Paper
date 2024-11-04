
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load the utility functions
source("Benchmarking/Simulation/utilities.R")

theme_set(theme_bw())

## GLOBAL VARIABLES ----
SETTING = 1
SAVE = FALSE
SHOW = TRUE
RUN = FALSE

# Import the sgdGMF package (only if you want to run the algorithms)
if (RUN) devtools::load_all()

# choose the method to be shown in the plots
.PEARSON = FALSE
.DEVIANCE = FALSE
.GLMPCA = TRUE
.AVAGRAD = TRUE
.FISHER = TRUE
.NBWAVE = TRUE
.CMF = TRUE
.NMF = TRUE
.NNLM = TRUE
.AIRWLS = TRUE
.NEWTON = TRUE
.SGD = TRUE

# set the path where to find the data
FILEPATH = paste("Benchmarking", "Simulation", "data", sep = "/")
IMGPATH = paste("Benchmarking", "Simulation", "img", sep = "/")

# set the custom colors for the clusters
myggcol = c("#F8766D", "#619CFF", "#00BA38", "#ff9933", "#d97ff2")
seaborn = c("#ff7c04", "#387cbc", "#e81c1c", "#50ac4c", "#a04ca4")

## PLOTTING FUNCTIONS ----
plot.tsne.grid = function (tsne, by = 1) {
  models = c()
  df = data.frame(
    x = c(),
    y = c(),
    group = c(),
    batch = c(),
    model = c())

  for (t in tsne) {
    n = nrow(t$tsne)
    dft = data.frame(
      x = t$tsne[,1],
      y = t$tsne[,2],
      group = t$group,
      batch = t$batch,
      model = rep(t$model, n))
    idx = seq(1, n, by = by)
    df = rbind(df, dft[idx, ])
    models = c(models, t$model)
  }

  df$x = as.numeric(df$x)
  df$y = as.numeric(df$y)
  df$group = as.factor(df$group)
  df$batch = as.factor(df$batch)
  df$model = factor(df$model, levels = models)

  levels(df$group) = 1:length(unique(df$group))
  levels(df$batch) = 1:length(unique(df$batch))

  gg_color_hue = function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group, pch = batch)) +
    geom_point(alpha = 0.8, size = 1.5) + facet_wrap(vars(model)) +
    labs(color = "Celltype", pch = "Batch") +
    scale_color_brewer(palette = "Set2") +
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 13),
          strip.text = element_text(size = 13),
          strip.background = element_rect(color = "white"),
          panel.grid = element_blank(),
          axis.title = element_blank())

#    theme(text = element_text(size = 20),
#          legend.position = "bottom",
#          legend.text = element_text(size = 10),
#          legend.title = element_text(size = 15),
#          axis.title = element_blank(),
#          axis.text = element_blank())

  return (plt)
}

plot.sil.score = function (sil, main = "Silhouette plot", by = 1) {
  n = nrow(sil)
  sil = cluster::sortSilhouette(sil)
  df = data.frame(index = 1:5000, cluster = as.factor(sil[,1]), width = sil[,3])
  mn = mean(df$width)
  df = df[seq(1, n, by = by),]
  plt = ggplot(data = df, map = aes(x = index, y = width, fill = cluster, color = cluster)) +
    geom_area() + ylim(-1,+1) +
    geom_hline(yintercept = mn, col = 2, lty = 2) +
    labs(x = "Cell index", y = "Silhouette width") +
    labs(color = "Cell-type", fill = "Cell-type") +
    labs(title = paste(main, " (", round(mn, 3), ") ", sep = ""))
  return(plt)
}

plot.sil.grid = function (sil, by = 1) {

  models = c()
  means = c()

  df = data.frame(
    model = c(),
    index = c(),
    cluster = c(),
    width = c(),
    mean = c())

  for (l in sil) {
    n = nrow(l$sil)
    s = cluster::sortSilhouette(l$sil)
    m = mean(s[,3])
    dft = data.frame(
      model = rep(l$model, n),
      index = 1:n,
      cluster = s[,1],
      width = s[,3],
      mean = rep(m, n))
    idx = seq(1, n, by = by)
    df = rbind(df, dft[idx, ])
    models = c(models, l$model)
    means = c(means, round(m, 3))
  }

  df$model = factor(df$model, levels = models)
  df$index = as.numeric(df$index)
  df$cluster = as.factor(df$cluster)
  df$width = as.numeric(df$width)
  df$mean = as.numeric(df$mean)

  txt = as.character(round(unique(df$mean), 3))
  plt = ggplot(data = df, map = aes(x = index, fill = cluster, color = cluster)) +
    geom_area(map = aes(y = width)) + geom_line(map = aes(y = mean), col = 1, lty = 2) +
    facet_wrap(vars(model), scales = "fixed") + ylim(-1,+1) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Cell index", y = "Silhouette width") +
    labs(color = "Cell-type", fill = "Cell-type") +
    theme(strip.background = element_blank())


  plt = plt + geom_text(
    data = data.frame(
      model = factor(models, levels = models),
      cluster = as.factor(rep(1, length(sil))),
      x = rep(10, length(sil)),
      y = means+0.15),
    map = aes(x = x, y = y),
    label = means,
    col = 1,
    hjust = .0,
    vjust = .5,
    size = 3
  )

  return(plt)
}

load.data = function (n, m, d, i) {

  setting = "bubble"
  filepath = FILEPATH
  fileid = paste("_s", setting, "_n", n, "_m", m, "_d", d, "_i", i, sep = "")
  filename = paste("summary", fileid, ".csv", sep = "")

  df = read.table(file = paste(filepath, filename, sep = "/"),
                  header = TRUE, dec = ".", sep = ";")

  return (df)
}

filter.data = function (df, n, m, d, i) {

  # Filter out unused rows
  df = df[-which(df$Set != "Test"), ]
  df = df[-which(df$Model == "C-SGD"), ]

  # Create new variables
  df = df %>% mutate(Dim = n/100, NComp = d)

  # Rename and sort the variables
  df = data.frame(
    iteration = df$Iter,
    model = df$Model,
    dimension = df$Dim,
    ncomp = df$NComp,
    time = df$Time,
    error = df$RSS,
    deviance = df$Dev,
    silhouette = df$Sil,
    purity = df$Purity)

  # Rename some model labels
  df$model[df$model == "NNLM"] = "NMF+"
  df$model[df$model == "B-SGD"] = "SGD"

  # Transform model into a factor variable
  models = c("CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "AIRWLS", "Newton", "SGD")
  df$model = factor(df$model, levels = models)

  return(df)
}

plot.summary = function (s = "bubble", n = 5000, m = 500, d = 5, i = 25) {

  require(ggh4x)
  require(ggbreak)

  df = load.data(n, m, d, i)
  df = filter.data(df, n, m, d, i)
  df = df[, -which(colnames(df) %in% c("iteration", "dimension", "ncomp"))]
  df$deviance = 100 * df$deviance
  df$error = 100 * df$error
  df = df[, -which(colnames(df) == "error")]

  df2 = reshape2::melt(df, id.var = "model")

  df2$variable = factor(df2$variable, levels = c("time", "deviance", "silhouette", "purity"))
  levels(df2$variable) = c("Time", "Deviance", "Silhouette", "Purity")

  scales <- list(
    scale_y_continuous(limits = range(df$time[df$model != "AvaGrad"], na.rm = TRUE)),
    scale_y_continuous(limits = range(df$deviance)),
    scale_y_continuous(limits = range(df$silhouette)),
    scale_y_continuous(limits = range(df$purity))
  )

  title = join.string("n = ", n, ", m = ", m, ", d = ", d)
  palette = "Temps"
  colors = hcl.colors(9, palette = "Temps")
  plt = ggplot(data = df2, map = aes(x = model, y = value, color = model)) +
    geom_boxplot(map = aes(fill = model), alpha = 0.5) +
    geom_boxplot(alpha = 0.5) +
    facet_grid(rows = vars(variable), scales = "free_y") +
    facetted_pos_scales(y = scales) + ggtitle(title) +
    # scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
    # scale_color_brewer(palette = palette) + scale_fill_brewer(palette = palette) +
    theme(axis.title = element_blank(), axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
          strip.text = element_text(size = 13), strip.background = element_rect(color = "white"), # element_blank(),
          plot.title = element_text(size = 13), legend.position = "none",
          legend.title = element_text(size = 13), legend.text = element_text(size = 10))

  return (plt)
}

## LOAD DATA ----
filename = "example_sbubble_n5000_m500_d5_i25.RData"
load(file = paste(FILEPATH, filename, sep = "/"))

## DATA EXTRACTION ----
logcounts = as.data.frame(logcounts(sim))
counts = as.data.frame(counts(sim))
cells = as.data.frame(colData(sim))
genes = as.data.frame(rowData(sim))
meta = metadata(sim)
groups = as.numeric(as.factor(cells$Group))
batches = as.numeric(as.factor(cells$Batch))

# PCA embedding
pca = RSpectra::svds(scale(t(as.matrix(logcounts))), k = 10)$u

# t-SNE embedding
tsne = Rtsne::Rtsne(pca, dims = 2, verbose = TRUE, num_threads = 8)$Y

df = data.frame(
  x = c(minmax(pca[,1]), minmax(tsne[,1])),
  y = c(minmax(pca[,2]), minmax(tsne[,2])),
  group = as.factor(rep(groups, 2)),
  batch = as.factor(rep(batches, 2)),
  embedding = as.factor(rep(c("PCA", "tSNE"), each = n))
)

ggplot(data = df, map = aes(x = x, y = y, color = group, pch = batch)) +
  geom_point(alpha = 0.99) + facet_grid(cols = vars(embedding)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "PC1", y = "PC2", color = "Cell-type", pch = "Batch") +
  theme(strip.background = element_blank())

## TRAIN-TEST SPLIT ----
X = model.matrix(~ Batch, data = cells)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)

test = Y
train = Y
ctrain = Y

test[data$test] = NA
train[data$train] = NA
ctrain = naive.completion(train)

## MODEL FIT ----

if (RUN) {
  model.pearson = fit.pearson(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family)
  model.deviance = fit.deviance(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family)
  model.avagrad = fit.avagrad(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 1000, tol = 1e-04)
  model.fisher = fit.fisher(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-05)
  model.nbwave = fit.nbwave(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-04)
  model.cmf = fit.cmf(y = train, x = X, z = NULL, ncomp = ncomp, family = family, verbose = FALSE, maxiter = 500)
  model.nmf = fit.nmf(y = ctrain, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE)
  model.nnlm = fit.nnlm(y = train, x = NULL, z = NULL, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000)
  model.airwls = fit.C.airwls(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.9)
  model.newton = fit.C.newton(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, stepsize = 0.2)
  model.sgd = fit.C.bsgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000, stepsize = 0.01)

  filename = "models_sbubble_n5000_m500_d5_i25.RData"
  save(model.pearson, model.deviance, model.glmpca, model.avagrad,
       model.fisher, model.nbwave, model.cmf, model.nmf, model.nnlm,
       model.airwls, model.newton, model.sgd,
       file = paste(FILEPATH, filename, sep = "/"))
} else {
  filename = "models_sbubble_n5000_m500_d5_i25.RData"
  load(file = paste(FILEPATH, filename, sep = "/"))
}



## TSNE PROJECTION ----
list.tsne = list()
if (.CMF) list.tsne$cmf = list(model = "CMF (cmfrec)", tsne = model.cmf$tsne, group = groups, batch = batches)
if (.NMF) list.tsne$nmf = list(model = "NMF (NMF)", tsne = model.nmf$tsne, group = groups, batch = batches)
if (.NNLM) list.tsne$nnlm = list(model = "NMF+ (NNLM)", tsne = model.nnlm$tsne, group = groups, batch = batches)
if (.PEARSON) list.tsne$pearson = list(model = "Pearson", tsne = model.pearson$tsne, group = groups, batch = batches)
if (.DEVIANCE) list.tsne$deviance = list(model = "Deviance", tsne = model.deviance$tsne, group = groups, batch = batches)
if (.AVAGRAD) list.tsne$avagrad = list(model = "AvaGrad (glmPCA)", tsne = model.avagrad$tsne, group = groups, batch = batches)
if (.FISHER) list.tsne$fisher = list(model = "Fisher (glmPCA)", tsne = model.fisher$tsne, group = groups, batch = batches)
if (.NBWAVE) list.tsne$nbwave = list(model = "NBWaVE (NewWave)", tsne = model.nbwave$tsne, group = groups, batch = batches)
if (.AIRWLS) list.tsne$airwls = list(model = "AIRWLS (sgdGMF)", tsne = model.airwls$tsne, group = groups, batch = batches)
if (.NEWTON) list.tsne$newton = list(model = "Newton (sgdGMF)", tsne = model.newton$tsne, group = groups, batch = batches)
if (.SGD) list.tsne$sgd = list(model = "SGD (sgdGMF)", tsne = model.sgd$tsne, group = groups, batch = batches)

plt.tsne = plot.tsne.grid(list.tsne, by = 5)

if (SHOW) print(plt.tsne)

if (SAVE) {
  filename = paste("simulation_example_tsne.pdf", sep = "")
  path = IMGPATH; width = 1.4; height = 1.6; zoom = 14
  ggsave(filename = filename, plot = plt.tsne, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}

## SILHOUETTE SCORES ----
list.sil = list()
if (.CMF) list.sil$cmf = list(model = "CMF (cmfrec)", sil = cluster::silhouette(groups, dist(model.cmf$tsne)))
if (.NMF) list.sil$nmf = list(model = "NMF (NMF)", sil = cluster::silhouette(groups, dist(model.nmf$tsne)))
if (.NNLM) list.sil$nnlm = list(model = "NMF+ (NNLM)", sil = cluster::silhouette(groups, dist(model.nnlm$tsne)))
if (.PEARSON) list.sil$pearson = list(model = "Pearson", sil = cluster::silhouette(groups, dist(model.pearson$tsne)))
if (.DEVIANCE) list.sil$deviance = list(model = "Deviance", sil = cluster::silhouette(groups, dist(model.deviance$tsne)))
if (.GLMPCA) list.sil$avagrad = list(model = "AvaGrad (glmPCA)", sil = cluster::silhouette(groups, dist(model.avagrad$tsne)))
if (.GLMPCA) list.sil$fisher = list(model = "Fisher (glmPCA)", sil = cluster::silhouette(groups, dist(model.fisher$tsne)))
if (.NBWAVE) list.sil$nbwave = list(model = "NBWaVE (NewWave)", sil = cluster::silhouette(groups, dist(model.nbwave$tsne)))
if (.AIRWLS) list.sil$airwls = list(model = "AIRWLS (sgdGMF)", sil = cluster::silhouette(groups, dist(model.airwls$tsne)))
if (.NEWTON) list.sil$newton = list(model = "Newton (sgdGMF)", sil = cluster::silhouette(groups, dist(model.newton$tsne)))
if (.SGD) list.sil$sgd = list(model = "SGD (sgdGMF)", sil = cluster::silhouette(groups, dist(model.sgd$tsne)))

plt.sil = plot.sil.grid(list.sil, by = 5)

if (SHOW) print(plt.sil)

if (SAVE) {
  filename = paste("simulation_example_silhouette.pdf", sep = "")
  path = IMGPATH; width = 1.4; height = 1.6; zoom = 14;
  ggsave(filename = filename, plot = plt.sil, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}

## LOAD SUMMARY STATS ----
plt.stat = plot.summary(n = 5000, m = 500, d = 5, i = 100)

plt = ggpubr::ggarrange(
  plt.stat + ggtitle("Summary statistics"),
  plt.tsne + ggtitle(" tSNE projections"),
  nrow = 1, ncol = 2, widths = c(1,2), labels = "AUTO")

if (SHOW) print(plt)

if (SAVE) {
  filename = paste("simulation_example.pdf", sep = "")
  path = IMGPATH; width = 2; height = 1.6; zoom = 14
  ggsave(filename = filename, plot = plt, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}

## END OF FILE ----
