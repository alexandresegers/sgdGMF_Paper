
## WORKSPACE SETUP ----

## Clean the workspace
rm(list = ls())
graphics.off()

## Load the package
devtools::load_all()
devtools::document()

## Load the utility functions
source("sim/utilities.R")

## GLOBAL VARIABLES ----
SETTING = 1
SAVE = TRUE
SHOW = TRUE

theme_set(theme_bw())

PALETTE = "Set2"

colors = c("#F8766D", "#619CFF", "#00BA38", "#ff9933", "#d97ff2")

## PLOTTING FUNCTIONS ----
summary.cv.gof = function (df) {

  require(dplyr)

  cv.med = df %>% dplyr::group_by(ncomp) %>%
    dplyr::summarise(AIC = median(aic), BIC = median(bic), Deviance = median(dev)) %>%
    reshape2::melt(id.vars = "ncomp") %>% as.data.frame()

  cv.min = df %>% dplyr::group_by(ncomp) %>%
    dplyr::summarise(AIC = min(aic), BIC = min(bic), Deviance = min(dev)) %>%
    reshape2::melt(id.vars = "ncomp") %>% as.data.frame()

  cv.max = df %>% dplyr::group_by(ncomp) %>%
    dplyr::summarise(AIC = max(aic), BIC = max(bic), Deviance = max(dev)) %>%
    reshape2::melt(id.vars = "ncomp") %>% as.data.frame()

  cv.gof = data.frame(
    ncomp = as.factor(cv.med$ncomp), variable = as.factor(cv.med$variable),
    value = cv.med$value, lower = cv.min$value, upper = cv.max$value)

  return (cv.gof)
}

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

  # colors = c("#ff7c04", "#387cbc", "#e81c1c", "#50ac4c", "#a04ca4")
  plt = ggplot(data = df, mapping = aes(x = x, y = y, color = group)) +
    geom_point(alpha = 0.5, size = 1.5) +
    facet_wrap(vars(model), ncol = 2, nrow = 3) +
    labs(color = "Cell-type", pch = "Batch") +
    scale_color_brewer(palette = PALETTE) +
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
    scale_color_brewer(palette = PALETTE) +
    scale_fill_brewer(palette = PALETTE) +
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
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
      model = rep(l$model, times = n),
      index = seq(n),
      cluster = as.numeric(s[,1]),
      width = as.numeric(s[,3]),
      mean = rep(m, times = n))
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
    scale_color_brewer(palette = PALETTE) +
    scale_fill_brewer(palette = PALETTE) +
    # scale_colour_manual(values = colors) +
    # scale_fill_manual(values = colors) +
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

plot.stat.summary = function (s = "bubble", n = 5000, m = 500, d = 5, i = 25) {

  fileid = join.string("_s", s, "_n", n, "_m", m, "_d", d, "_i", i)
  filepath = join.path("data", "splatter")
  filename = join.string("summary", fileid, ".csv")

  models = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SGD", "B-SGD")
  variables = c("Time", "Deviance", "Error", "Silhouette")

  df = read.table(file = join.path(filepath, filename), header = TRUE, sep = ";")
  df = df[df$Set == "Test", ]
  df = df[, c("Model", "Time", "Dev", "RSS", "Sil")]

  df = data.frame(
    Model = factor(rep(df$Model, times = 4), levels = models),
    Val = c(log10(df$Time), 100 * df$Dev, 100 * df$RSS, df$Sil),
    Var = factor(rep(variables, each = nrow(df)), levels = variables))

  levels(df$Model) = c("CMF", "NMF", "NNLM", "glmPCA", "NBWaVE", "AIRWLS", "Newton", "C-SQN", "B-SQN")

  plt = ggplot(data = df, map = aes(x = Model, y = Val, color = Model, fill = Model)) +
    geom_boxplot(alpha = 0.5) + facet_grid(rows = vars(Var), scales = "free_y") +
    ggtitle(join.string("n = ", 5000, ", m = ", 500, ", d = ", 5)) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
          strip.text = element_text(size = 13),
          plot.title = element_text(size = 13),
          legend.position = "none",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 10))

  return (plt)
}

## LOAD DATA ----
filename = "example_sbubble_n5000_m500_d5_i25.RData"
filepath = join.path("data", "splatter", filename)
load(file = filepath)

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
  geom_point(alpha = 0.9) + facet_grid(cols = vars(embedding)) +
  scale_color_brewer(palette = PALETTE) +
  labs(x = "PC1", y = "PC2", color = "Cell-type", pch = "Batch")

## TRAIN-TEST SPLIT ----
X = model.matrix(~ Batch, data = cells)
Z = matrix(1, nrow = m, ncol = 1)
Y = matrix(NA, nrow = n, ncol = m)
Y[] = t(counts)

## MODEL SELECTION ----
method = "bsgd"
ctr.cv = list(nfolds = 20, verbose = TRUE, parallel = FALSE, nthreads = 4)
ctr.alg = list(maxiter = 1000, size = c(100,20), frequency = 250, rate0 = 0.01)

RUN = FALSE
if (RUN) {
  gmf.fit = sgdGMF::sgdgmf.cv(
    Y, X, Z, ncomp = 1:10, family = family, method = method,
    control.alg = ctr.alg, control.cv = ctr.cv)

  cv.gof = summary.cv.gof(gmf.fit$summary.cv)
} else {
  cv.gof = read.table(
    file = "data/splatter/example_sbubble_n5000_m500_d5_i25.csv",
    header = TRUE, sep = ";", dec = ".")

  cv.gof = summary.cv.gof(cv.gof)
}

plt.gof = cv.gof %>%
  ggplot(map = aes(x = as.factor(ncomp), y = value, ymin = lower, ymax = upper)) +
  geom_errorbar(color = "grey50", width = 0.2) +
  geom_point(color = "red", size = 2) +
  facet_grid(rows = vars(as.factor(variable)), scale = "free_y") +
  labs(x = "Latent space rank", y = "Goodness of fit") +
  ggtitle("Model selection criteria") +
  theme(plot.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        axis.title.y = element_blank())

print(plt.gof)

## SEQUENTIAL MODEL FIT ----
gmf.fit.list = list()
for (ncomp in 1:10) {
  cat(" Rank =", ncomp, "\n")
  # Fit the model
  fit = sgdGMF::sgdgmf.fit(
    Y, X, Z, ncomp = ncomp, family = family,
    method = method, control.alg = ctr.alg)

  # Store the estimated model
  gmf.fit.list[[ncomp]] = fit
}

for (ncomp in 1:10) {
  cat(" Rank =", ncomp, "\n")
  fit = gmf.fit.list[[ncomp]]

  # Compute a 2D tSNE embedding
  tsn = scale(Rtsne::Rtsne(as.matrix(fit$U), dims = 2,
                           parallel = TRUE, num_threads = 8)$Y)

  # Compute the silhouette obs-by-obs
  sil = cluster::silhouette(groups, dist(tsn))

  dfsil = as.data.frame(cbind(ncomp, sil[,-2]))
  colnames(dfsil) = c("ncomp", "group", "value")
  dfsil$ncomp = as.factor(dfsil$ncomp)
  dfsil$group = as.factor(dfsil$group)

  # Average the silhouette group-by-group
  avgsil = dfsil %>% group_by(ncomp, group) %>%
    reframe(value = mean(value)) %>% as.data.frame()

  # Compute the purity obs-by-obs
  u = prcomp(tcrossprod(fit$U, fit$V))$x[,1:ncomp]
  pur = bluster::neighborPurity(u, clusters = groups)

  dfpur = as.data.frame(cbind(ncomp, groups, pur$purity))
  colnames(dfpur) = c("ncomp", "group", "value")
  dfpur$ncomp = as.factor(dfpur$ncomp)
  dfpur$group = as.factor(dfpur$group)

  # Average the silhouette group-by-group
  avgpur = dfpur %>% group_by(ncomp, group) %>%
    reframe(value = mean(value)) %>% as.data.frame()

  # Store the results
  fit$tsne = tsn
  fit$sil = sil
  fit$pur = pur
  fit$avgsil = avgsil
  fit$avgpur = avgpur
  gmf.fit.list[[ncomp]] = fit
}


## SUMMARY PLOTS ----

# Silhouette summary
df.sil = do.call("rbind", lapply(gmf.fit.list, function(x) x$avgsil))
df.sil$variable = as.factor(rep("Silhouette", times = 50))
plt.sil = df.sil %>%
  ggplot(map = aes(x = as.numeric(ncomp), y = value, color = group)) +
  geom_line() + geom_point() + facet_grid(rows = vars(variable)) +
  labs(x = "Latent space rank", y = "Silhouette") +
  scale_color_brewer(palette = PALETTE) +
  scale_x_continuous(breaks = 1:10) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.text = element_text(size = 13),
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank())

print(plt.sil)

# Purity summary
df.pur = do.call("rbind", lapply(gmf.fit.list, function(x) x$avgpur))
df.pur$variable = as.factor(rep("Purity", times = 50))
plt.pur = df.pur %>%
  ggplot(map = aes(x = as.numeric(ncomp), y = value, color = group)) +
  geom_line() + geom_point() + facet_grid(rows = vars(variable)) +
  labs(x = "Latent space rank", y = "Purity") +
  scale_color_brewer(palette = PALETTE) +
  scale_x_continuous(breaks = 1:10) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.text = element_text(size = 13),
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank())

print(plt.pur)

# tSNE embeddings
plt.tsne = plot.tsne.grid(list(
  list(model = "Matrix rank: 2", tsne = gmf.fit.list[[2]]$tsne, group = groups, batch = batches),
  list(model = "Matrix rank: 3", tsne = gmf.fit.list[[3]]$tsne, group = groups, batch = batches),
  list(model = "Matrix rank: 4", tsne = gmf.fit.list[[4]]$tsne, group = groups, batch = batches),
  list(model = "Matrix rank: 5", tsne = gmf.fit.list[[5]]$tsne, group = groups, batch = batches),
  list(model = "Matrix rank: 6", tsne = gmf.fit.list[[6]]$tsne, group = groups, batch = batches),
  list(model = "Matrix rank: 10", tsne = gmf.fit.list[[10]]$tsne, group = groups, batch = batches)
), by = 5)

print(plt.tsne)

# Silhouette profile
plt.sil = plot.sil.grid(list(
  list(model = "rank = 2", sil = gmf.fit.list[[2]]$sil),
  list(model = "rank = 3", sil = gmf.fit.list[[3]]$sil),
  list(model = "rank = 4", sil = gmf.fit.list[[4]]$sil),
  list(model = "rank = 5", sil = gmf.fit.list[[5]]$sil),
  list(model = "rank = 6", sil = gmf.fit.list[[6]]$sil),
  list(model = "rank = 10", sil = gmf.fit.list[[10]]$sil)
), by = 5)

print(plt.sil)

# Summary
plt.sil = rbind(df.sil, df.pur) %>%
  ggplot(map = aes(x = as.numeric(ncomp), y = value, color = group)) +
  geom_line() + geom_point() + facet_grid(rows = vars(variable), scale = "free") +
  labs(x = "Latent space rank", y = "") +
  scale_color_brewer(palette = PALETTE) +
  scale_x_continuous(breaks = 1:10) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.text = element_text(size = 13),
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank())

plt.gof = cv.gof %>%
  ggplot(map = aes(x = as.factor(ncomp), y = value, ymin = lower, ymax = upper)) +
  geom_errorbar(color = "grey50", width = 0.2) +
  geom_point(color = "red", size = 2) +
  facet_grid(rows = vars(as.factor(variable)), scale = "free_y") +
  labs(x = "Latent space rank", y = "Goodness of fit") +
  ggtitle("Goodness of fit criteria") +
  theme(plot.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

plt.stat = ggpubr::ggarrange(
  plt.gof, plt.sil, nrow = 2, heights = c(2.5,2), align = "v")

plt = ggpubr::ggarrange(
  plt.stat + ggtitle("Model selection criteria"),
  plt.tsne + ggtitle("tSNE projection"),
  nrow = 1, ncol = 2, widths = c(2,3))

print(plt)

if (SAVE) {
  filename = paste("model_selection.pdf", sep = "")
  path = "img/splatter"
  width = 2
  # height = 1.25
  height = 1.75
  zoom = 13
  ggsave(filename = filename, plot = plt, path = path,
         width = zoom * width, height = zoom * height, units = "cm")
}
