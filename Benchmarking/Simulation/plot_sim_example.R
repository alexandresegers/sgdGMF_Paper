
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load the package
devtools::load_all()

## Load the utility functions
source("sim/utilities.R")

theme_set(theme_bw())

## GLOBAL VARIABLES ----
SETTING = 1
SAVE = TRUE
SHOW = TRUE

colors = c("#F8766D", "#619CFF", "#00BA38", "#ff9933", "#d97ff2")

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

  colors = c("#ff7c04", "#387cbc", "#e81c1c", "#50ac4c", "#a04ca4")
  plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group, pch = batch)) +
    geom_point(alpha = 0.8, size = 2) + facet_wrap(vars(model)) +
    labs(color = "Cell-type", pch = "Batch") +
    scale_color_brewer(palette = "Set2") +
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 13),
          strip.text = element_text(size = 13),
          axis.title = element_blank())

#    theme(text = element_text(size = 20),
#          legend.position = "bottom",
#          legend.text = element_text(size = 10),
#          legend.title = element_text(size = 15),
#          axis.title = element_blank(),
#          axis.text = element_blank())

  if (SAVE) {
    filename = paste("example_sbubble_n5000_m500_d5_i25.pdf", sep = "")
    path = "img/splatter"
    width = 1
    height = 1
    zoom = 8
    ggsave(filename = filename, plot = plt, path = path,
           width = zoom * width, height = zoom * height, units = "cm")
  }

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
    labs(color = "Cell-type", fill = "Cell-type")


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
  file.path = paste("data", "bubble", sep = "/")
  file.id = paste("_s", setting, "_n", n, "_m", m, "_d", d, "_i", i, sep = "")
  file.name = paste("summary", file.id, ".csv", sep = "")

  df = read.table(file = paste(file.path, file.name, sep = "/"),
                  header = TRUE, dec = ".", sep = ";")

  models = c("CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "AIRWLS", "Newton", "SGD")
  df$model = factor(df$model, levels = models)

  return (df)
}

plot.stat.summary = function (s = "bubble", n = 5000, m = 500, d = 5, i = 25) {

  require(ggh4x)
  require(ggbreak)

  df = load.data(n = 5000, m = 500, d = 5, i = 25)
  df = df[, -which(colnames(df) %in% c("iteration", "dimension", "ncomp"))]
  df$deviance = 100 * df$deviance
  df$error = 100 * df$error

  df2 = reshape2::melt(df, id.var = "model")

  df2$variable = factor(df2$variable, levels = c("time", "deviance", "error", "silhouette"))
  levels(df2$variable) = c("Time", "Deviance", "Error", "Silhouette")

  scales <- list(
    scale_y_continuous(limits = c(0, 125)),
    scale_y_continuous(limits = range(df$deviance)),
    scale_y_continuous(limits = range(df$error)),
    scale_y_continuous(limits = range(df$silhouette))
  )

  title = join.string("n = ", 5000, ", m = ", 500, ", d = ", 5)

  plt = ggplot(data = df2, map = aes(x = model, y = value, color = model, fill = model)) +
    geom_boxplot(alpha = 0.5) + facet_grid(rows = vars(variable), scales = "free_y") +
    facetted_pos_scales(y = scales) + ggtitle(title) +
    theme(axis.title = element_blank(), axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
          strip.text = element_text(size = 13), plot.title = element_text(size = 13),
          legend.title = element_text(size = 13), legend.text = element_text(size = 10),
          legend.position = "none")

  return (plt)
}

## LOAD DATA ----
file.name = "example_sbubble_n5000_m500_d5_i25.RData"
file.path = join.path("data", "splatter", file.name)
load(file = file.path)

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
  labs(x = "PC1", y = "PC2", color = "Cell-type", pch = "Batch")

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
model.pearson = fit.pearson(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family)
model.glmpca = fit.glmpca(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-05)
model.nbwave = fit.nbwave(y = ctrain, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 200, tol = 1e-04)
model.sgd = fit.C.bsgd(y = train, x = X, z = Z, ncomp = ncomp, family = family, verbose = TRUE, maxiter = 2000, stepsize = 0.01)

## TSNE PROJECTION ----
plt.tsne = plot.tsne.grid(list(
  list(model = "Pearson", tsne = model.pearson$tsne, group = groups, batch = batches),
  list(model = "glmPCA", tsne = model.glmpca$tsne, group = groups, batch = batches),
  list(model = "NBWaVE", tsne = model.nbwave$tsne, group = groups, batch = batches),
  list(model = "SGD", tsne = model.sgd$tsne, group = groups, batch = batches)
), by = 5)

if (SHOW) {
  print(plt.tsne)
}

## SILHOUETTE SCORES ----
sil.pearson = cluster::silhouette(groups, dist(model.pearson$tsne))
sil.glmpca = cluster::silhouette(groups, dist(model.glmpca$tsne))
sil.nbwave = cluster::silhouette(groups, dist(model.nbwave$tsne))
sil.sgd = cluster::silhouette(groups, dist(model.sgd$tsne))

plt.sil = plot.sil.grid(list(
  list(model = "Pearson", sil = sil.pearson),
  list(model = "glmPCA", sil = sil.glmpca),
  list(model = "NBWaVE", sil = sil.nbwave),
  list(model = "SGD", sil = sil.sgd)
), by = 5)

if (SHOW) {
  print(plt.sil)
}

## LOAD SUMMARY STATS ----
df = load.data(n = 5000, m = 500, d = 5, i = 25)

plt.stat = plot.stat.summary()

plt = ggpubr::ggarrange(
  plt.stat + ggtitle("Summary statistics"),
  plt.tsne + ggtitle("tSNE projection"),
  nrow = 1, ncol = 2, widths = c(1,2))

if (SHOW) {
  print(plt)
}

if (SAVE) {
  filename = paste("example_sbubble_n5000_m500_d5_i25.pdf", sep = "")
  path = "img/splatter"
  width = 2
  height = 1.25
  zoom = 14
  ggsave(filename = filename, plot = plt, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}

## END OF FILE ----
