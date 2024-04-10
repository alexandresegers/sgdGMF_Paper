
## WORKSPACE SETUP ----

# Clean the workspace
rm(list = ls())
graphics.off()
options(warn = -1)

# Packages
library(ggplot2)
library(dplyr)
library(latex2exp)

theme_set(theme_bw())

# Global variables
SAVE = TRUE
SHOW = TRUE

## PLOTTING FUNCTIONS ----

load.data.old = function (n, m, d, i, write = FALSE) {

  file.path = paste("data", "splatter", sep = "/")
  file.id = paste("_sbubble", "_n", n, "_m", m, "_d", d, "_i", i, sep = "")
  file.old = paste("summary", file.id, ".csv", sep = "")
  file.new = paste("summary_glmpca", file.id, ".csv", sep = "")

  df.old = read.table(file = paste(file.path, file.old, sep = "/"), header = TRUE, dec = ".", sep = ";")
  df.new = read.table(file = paste(file.path, file.new, sep = "/"), header = TRUE, dec = ".", sep = ";")

  df = rbind(df.old, df.new)
  df = df[ which(df$Set   ==   "Test"), ]
  df = df[-which(df$Model == "glmPCA"), ]
  df = df[-which(df$Model ==  "C-SGD"), ]
  df = df[, -c(4, 7)]
  colnames(df) = c("model", "time", "error", "deviance", "silhouette", "iteration")

  models = c("CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "AIRWLS", "Newton", "SGD")
  df$model[df$model == "NNLM"] = "NMF+"
  df$model[df$model == "B-SGD"] = "SGD"
  df$model = factor(df$model, levels = models)
  df$dimension = floor(m / 10)
  df$ncomp = floor(d)

  df$time[is.infinite(df$time)] = NA
  df$error[is.infinite(df$error)] = NA
  df$deviance[is.infinite(df$deviance)] = NA

  if (write) {
    df %>%
      group_by(iteration, model, dimension, ncomp) %>%
      summarise(time, error, deviance, silhouette) %>%
      as.data.frame() %>%
      write.table(file = paste("data", "bubble", file.old, sep = "/"),
                  sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)
  }

  return (df)
}

load.data = function (n, m, d, i, write = FALSE) {

  # Define the file path and name
  file.path = paste("sim", "data", sep = "/")
  file.id = paste("_sbubble", "_n", n, "_m", m, "_d", d, "_i", i, sep = "")
  file.name = paste("summary", file.id, ".csv", sep = "")

  # Load the data-set
  df = read.table(file = paste(file.path, file.name, sep = "/"),
                  header = TRUE, dec = ".", sep = ";")

  models = c("CMF", "NMF", "NMF+", "AvaGrad",
             "Fisher", "NBWaVE", "AIRWLS", "Newton", "SGD")

  df$model = factor(df$model, levels = models)
  # df$dimension = floor(m / 10)
  # df$ncomp = floor(d)

  df$time[is.infinite(df$time)] = NA
  df$error[is.infinite(df$error)] = NA
  df$deviance[is.infinite(df$deviance)] = NA

  return (df)
}

group.data = function (df, melt = FALSE) {
  df = df %>%
    group_by(model, dimension, ncomp) %>%
    summarise(time = mean(time, na.rm = TRUE),
              error = mean(100 * error, na.rm = TRUE),
              deviance = mean(100 * deviance, na.rm = TRUE),
              silhouette = mean(silhouette, na.rm = TRUE)) %>%
    as.data.frame()

  if (melt){
    id = c("model", "dimension", "ncomp")
    df = df[, -which(colnames(df) == "error")]
    df = reshape2::melt(df, id.vars = id)
    levels(df$variable) = c("Time", "Deviance", "Silhouette")
  }

  return (df)
}

merge.data = function (df1, df2, melt = FALSE) {

  # Introduce an indicator for the simulation setting
  df1$setting = factor("A", levels = c("A", "B"))
  df2$setting = factor("B", levels = c("A", "B"))

  # Merge the two data-frames
  df = rbind(df1, df2)[, -6]

  # Group the data by model, matrix dimension and latent space dimension
  df = df %>%
    group_by(model, setting, dimension, ncomp) %>%
    summarise(time = mean(time, na.rm = TRUE),
              error = mean(error, na.rm = TRUE),
              deviance = mean(deviance, na.rm = TRUE),
              silhouette = mean(silhouette, na.rm = TRUE)) %>%
    as.data.frame()

  # Do you want to melt the dataset
  if (melt) {
    idxA = which(df$setting == "A")
    idxB = which(df$setting == "B")
    vars = c("Time", "Deviance", "Silhouette")
    dims = rep(NA, times = nrow(df))
    dims[idxA] = df$dimension[idxA]
    dims[idxB] = df$ncomp[idxB]

    df$dimension = dims
    df = df[, -which(colnames(df) == "ncomp")]
    df = df[, -which(colnames(df) == "error")]

    df = reshape2::melt(df, id.vars = c("model", "setting", "dimension"))
  }

  return (df)
}

melt.data = function (df) {

  idxA = which(df$setting == "A")
  idxB = which(df$setting == "B")
  vars = c("Time", "Deviance", "Silhouette")
  dims = rep(NA, times = nrow(df))
  dims[idxA] = df$dimension[idxA]
  dims[idxB] = df$ncomp[idxB]

  df$dimension = dims
  # df$time = log10(df$time)
  df = df[, -which(colnames(df) == "ncomp")]
  df = df[, -which(colnames(df) == "error")]

  df.melt = reshape2::melt(df, id.vars = c("model", "setting", "dimension"))

  return (df.melt)
}

## DATA LOAD AND PLOTS ----

# Load data
df.flat = data.frame()
for (r in c(10, 25, 50, 75, 100)) {
  df.temp = load.data(n = 100*r, m = 10*r, d = 5, i = 25, write = TRUE)
  df.flat = rbind(df.flat, df.temp)
}

df.deep = data.frame()
for (d in c(5, 10, 15, 20, 25)) {
  df.temp = load.data(n = 5000, m = 500, d = d, i = 25, write = TRUE)
  df.deep = rbind(df.deep, df.temp)
}

# Group the data by setting, model and dimensions
df.flat.group = group.data(df.flat, melt = TRUE)
df.deep.group = group.data(df.deep, melt = TRUE)

# Define the plot theme options
my.thm = theme(
  legend.position = "right", legend.title = element_text(size = 13),
  legend.text = element_text(size = 10), plot.title = element_text(size = 13),
  strip.text = element_text(size = 13), axis.title.y = element_blank())

# Plot the summary statistics under setting A
plt.flat = df.flat.group %>%
  ggplot(map = aes(x = dimension, y = value, color = model, fill = model)) +
  geom_line() + geom_point(size = 3) +
  labs(x = "Data matrix dimension", color = "Model", fill = "Model") +
  ggtitle(TeX("Setting A ($n \\uparrow, m \\uparrow, d = 5$)")) +
  facet_grid(rows = vars(variable), scale = "free_y") + my.thm

# Plot the summary statistics under setting B
plt.deep = df.deep.group %>%
  ggplot(map = aes(x = ncomp, y = value, color = model, fill = model)) +
  geom_line() + geom_point(size = 3) +
  labs(x = "Latent space rank", color = "Model", fill = "Model") +
  ggtitle(TeX("Setting B ($n = 5000, m = 500, d \\uparrow$)")) +
  facet_grid(rows = vars(variable), scale = "free_y") + my.thm

# Stack the two plots by column
plt = ggpubr::ggarrange(plt.flat, plt.deep, nrow = 1, ncol = 2,
                        common.legend = TRUE, legend = "right")

# Show the plot
if (SHOW) print(plt)

# Save the result
if (SAVE) {
  setting = "bubble"
  filepath = paste("img", "splatter", sep = "/")
  filename = paste("summary_", setting, "_sim.pdf", sep = "")
  zoom = 4
  nrows = 3
  ncols = 2
  width = 1 * ncols
  height = .6 * nrows
  ggplot2::ggsave(
    file = filename, path = filepath, plot = plt,
    width = zoom * width, height = zoom * height)
}

## ALTERNATIVE PLOT ----

# Merge the two simulation settings
df.melt = merge.data(df.flat, df.deep, melt = TRUE)

# This is an alternative plot with common scales
df.melt %>% ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point() +
  facet_grid(rows = vars(variable), cols = vars(setting), scale = "free") +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        axis.title = element_blank())

## INDIVIDUAL PLOTS ---

df.deep.group = df.deep %>%
  group_by(model, dimension, ncomp) %>%
  summarise(time = mean(time, na.rm = TRUE),
            error = mean(error, na.rm = TRUE),
            deviance = mean(deviance, na.rm = TRUE),
            silhouette = mean(silhouette, na.rm = TRUE))

df.flat.group = df.flat %>%
  group_by(model, dimension, ncomp) %>%
  summarise(time = mean(time, na.rm = TRUE),
            error = mean(error, na.rm = TRUE),
            deviance = mean(deviance, na.rm = TRUE),
            silhouette = mean(silhouette, na.rm = TRUE))

plt.deep.dev = df.deep.group %>% ggplot(map = aes(x = ncomp, y = deviance, color = model)) + geom_line() + geom_point(size = 3)
plt.deep.err = df.deep.group %>% ggplot(map = aes(x = ncomp, y = error, color = model)) + geom_line() + geom_point(size = 3)
plt.deep.sil = df.deep.group %>% ggplot(map = aes(x = ncomp, y = silhouette, color = model)) + geom_line() + geom_point(size = 3)
plt.deep.time = df.deep.group %>% ggplot(map = aes(x = ncomp, y = time, color = model)) + geom_line() + geom_point(size = 3)
plt.flat.dev = df.flat.group %>% ggplot(map = aes(x = dimension, y = deviance, color = model)) + geom_line() + geom_point(size = 3)
plt.flat.err = df.flat.group %>% ggplot(map = aes(x = dimension, y = error, color = model)) + geom_line() + geom_point(size = 3)
plt.flat.sil = df.flat.group %>% ggplot(map = aes(x = dimension, y = silhouette, color = model)) + geom_line() + geom_point(size = 3)
plt.flat.time = df.flat.group %>% ggplot(map = aes(x = dimension, y = time, color = model)) + geom_line() + geom_point(size = 3)

ggpubr::ggarrange(
  plt.flat.time, plt.deep.time,
  plt.flat.dev, plt.deep.dev,
  plt.flat.sil, plt.deep.sil,
  nrow = 3, ncol = 2, common.legend = TRUE, legend = "bottom")

ggpubr::ggarrange(
  plt.flat.time, plt.flat.dev, plt.flat.sil,
  plt.deep.time, plt.deep.dev, plt.deep.sil,
  nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom")

## END OF FILE ----

