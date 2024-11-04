
## WORKSPACE SETUP ----

# Clean the workspace
rm(list = ls())
graphics.off()
options(warn = -1)

# Packages
library(ggplot2)
library(dplyr)
library(latex2exp)

# Default theme
theme_set(theme_bw())

# Global variables
SAVE = FALSE
SHOW = TRUE

# File and image paths
FILEPATH = paste("Benchmarking", "Simulation", "data", sep = "/")
IMGPATH = paste("Benchmarking", "Simulation", "img", sep = "/")

# Color palette
COLORS = hcl(h = seq(15, 375, length = 10), l = 65, c = 100)[1:9]

## UTILITY FUNCTIONS ----

# Loads the data corresponding to:
# - n: data matrix rows
# - m: data matrix columns
# - d: fixed latent space rank
# - i: number of iterations in the simulation experiment
load.data = function (n, m, d, i) {

  # Define the file path and name
  filepath = FILEPATH
  fileid = paste("_sbubble", "_n", n, "_m", m, "_d", d, "_i", i, sep = "")
  filename = paste("summary", fileid, ".csv", sep = "")

  # Load the data-set
  df = read.table(file = paste(filepath, filename, sep = "/"),
                  header = TRUE, dec = ".", sep = ";")

  # Set the model dimensions
  df$Dim = m/10
  df$NComp = d

  # Return the transformed object
  return (df)
}

# Preprocesses the loaded data selecting, renaming and recasting
# the variables, excluding NA, NaN and infinite values
filter.data = function (df) {

  # define the new labels for the column names and for the model levels
  old.names = c("Model", "Time", "RSS", "Dev", "Sil", "Purity", "Iter", "Dim", "NComp")
  new.names = c("model", "time", "error", "deviance", "silhouette", "purity", "iteration", "dimension", "rank")
  old.models = c("CMF", "NMF", "NNLM", "AvaGrad", "Fisher", "NBWaVE", "AIRWLS", "Newton", "B-SGD")
  new.models = c("CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "AIRWLS", "Newton", "SGD")

  # Filter out the unused rows and columns
  # (do not change the order of the following lines of code)
  df = df[-which(df$Set != "Test"),]
  df = df[-which(df$Model == "C-SGD"),]
  df = df[,-which(colnames(df) == "Set")]
  df = df[,-which(colnames(df) == "Cos")]

  # Map the column names to the new labels
  colnames(df) <- new.names

  # Map the model levels to the new labels
  df$model = factor(df$model, levels = old.models)
  levels(df$model) = new.models

  # Fill the infinite values with NAs
  df$time[is.infinite(df$time)] = NA
  df$error[is.infinite(df$error)] = NA
  df$deviance[is.infinite(df$deviance)] = NA
  df$silhouette[is.infinite(df$silhouette)] = NA
  df$purity[is.infinite(df$purity)] = NA

  # Rebuild the data-frame changing the order of the columns
  df = data.frame(
    "iteration" = df$iteration,
    "model" = df$model,
    "dimension" = df$dimension,
    "rank" = df$rank,
    "time" = df$time,
    "error" = df$error,
    "deviance" = df$deviance,
    "silhouette" = df$silhouette,
    "purity" = df$purity)

  # Return the transformed data-frame
  return (df)
}

# Groups the data by model. matrix dimension and latent space rank
group.data = function (df, melt = FALSE) {
  df = df %>%
    group_by(model, dimension, rank) %>%
    summarise(time = mean(time, na.rm = TRUE),
              error = mean(100 * error, na.rm = TRUE),
              deviance = mean(100 * deviance, na.rm = TRUE),
              silhouette = mean(silhouette, na.rm = TRUE),
              purity = mean(purity, na.rm = TRUE)) %>%
    as.data.frame()

  if (melt){
    id = c("model", "dimension", "rank")
    df = df[, -which(colnames(df) == "error")]
    df = reshape2::melt(df, id.vars = id)
    levels(df$variable) = c("Time", "Deviance", "Silhouette", "Purity")
  }

  return (df)
}

# Merges and reshapes the data to obtain a melted data-frame
# which can fit in ggplot facet_grid o facet_wrap
merge.data = function (df1, df2, melt = FALSE) {

  # Introduce an indicator for the simulation setting
  df1$setting = factor("A", levels = c("A", "B"))
  df2$setting = factor("B", levels = c("A", "B"))

  # Merge the two data-frames
  df = rbind(df1, df2) # [, -6]

  # Group the data by model, matrix dimension and latent space dimension
  df = df %>%
    group_by(model, setting, dimension, rank) %>%
    summarise(time = mean(time, na.rm = TRUE),
              error = mean(error, na.rm = TRUE),
              deviance = mean(deviance, na.rm = TRUE),
              silhouette = mean(silhouette, na.rm = TRUE),
              purity = mean(purity, na.rm = TRUE)) %>%
    as.data.frame()

  # Do you want to melt the dataset
  if (melt) {
    idxA = which(df$setting == "A")
    idxB = which(df$setting == "B")
    vars = c("Time", "Deviance", "Silhouette", "Purity")
    dims = rep(NA, times = nrow(df))
    dims[idxA] = df$dimension[idxA]
    dims[idxB] = df$rank[idxB]

    df$dimension = dims
    df = df[, -which(colnames(df) == "rank")]
    df = df[, -which(colnames(df) == "error")]

    df = reshape2::melt(df, id.vars = c("model", "setting", "dimension"))
  }

  return (df)
}

melt.data = function (df) {

  idxA = which(df$setting == "A")
  idxB = which(df$setting == "B")
  vars = c("Time", "Deviance", "Silhouette", "Purity")
  dims = rep(NA, times = nrow(df))
  dims[idxA] = df$dimension[idxA]
  dims[idxB] = df$rank[idxB]

  df$dimension = dims
  # df$time = log10(df$time)
  df = df[, -which(colnames(df) == "rank")]
  df = df[, -which(colnames(df) == "error")]

  df.melt = reshape2::melt(df, id.vars = c("model", "setting", "dimension"))

  return (df.melt)
}

## DATA LOAD AND PREPROCESSING ----

# Flat data: fixed rank, increasing dimensions
df.flat = data.frame()
for (r in c(10, 25, 50, 75, 100)) {
  df.temp = load.data(n = 100*r, m = 10*r, d = 5, i = 100)
  df.temp = filter.data(df.temp)
  df.flat = rbind(df.flat, df.temp)
}

# Deep data: increasing rank, fixed dimensions
df.deep = data.frame()
for (d in c(5, 10, 15, 20, 25)) {
  df.temp = load.data(n = 5000, m = 500, d = d, i = 100)
  df.temp = filter.data(df.temp)
  df.deep = rbind(df.deep, df.temp)
}

# Set the plot theme options
my.thm = theme(
  legend.position = "right", legend.title = element_text(size = 13),
  legend.text = element_text(size = 10), plot.title = element_text(size = 13),
  strip.text = element_text(size = 13), axis.title.y = element_blank(),
  strip.background = element_rect(color = "white"))


## PLOT: FIRST LAYOUT ----

# Merge the two simulation settings
df.melt = merge.data(df.flat, df.deep, melt = TRUE)

# Map the variable labels to upper case letters
df.melt$variable = factor(df.melt$variable,
                          levels = c("time", "deviance", "silhouette", "purity"),
                          labels = c("Time", "Deviance", "Silhouette", "Purity"))

# Map the setting labels to TeX expressions
df.melt$setting = factor(df.melt$setting,
                         levels = c("A", "B"),
                         labels = c(TeX(paste("Setting A ($n \\uparrow, m \\uparrow, d = 5$)")),
                                    TeX(paste("Setting B ($n = 5000, m = 500, d \\uparrow$)"))))

# This is an alternative plot with common scales
plt = df.melt %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  facet_grid(rows = vars(variable), cols = vars(setting),
             scale = "free", labeller = label_parsed) +
  scale_colour_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(color = "Method") +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        strip.background = element_rect(color = "white"),
        strip.text = element_text(size = 13),
        axis.title = element_blank())

# Show the plot
if (SHOW) print(plt)

# Save the result
if (SAVE) {
  filepath = IMGPATH
  filename = "simulation_summary_4x2.pdf"
  zoom = 3.8
  nrows = 4
  ncols = 2
  width = 1 * ncols
  height = .6 * nrows
  ggplot2::ggsave(
    file = filename, path = filepath, plot = plt,
    width = zoom * width, height = zoom * height)
}

## PLOT: SECOND LAYOUT ----

# Separate the time summary to the other gof measures
df.time = df.melt[df.melt$variable == "Time", ]
df.gof = df.melt[df.melt$variable != "Time", ]

# Separate the execution times of the tow settings
df.A.time = df.melt[df.melt$variable == "Time" & df.melt$setting == levels(df.melt$setting)[1], ]
df.B.time = df.melt[df.melt$variable == "Time" & df.melt$setting == levels(df.melt$setting)[2], ]

# Execution time in the first setting
plt.A.time = df.A.time %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  scale_colour_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(x = "Matrix dimension", y = "Execution time (s)", color = "Method",
       title = TeX(paste("Setting A ($n \\uparrow, m \\uparrow, d = 5$)"))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

# Execution time in the second setting
plt.B.time = df.B.time %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  scale_colour_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(x = "Latent space rank", y = "Execution time (s)", color = "Method",
       title = TeX(paste("Setting B ($n = 5000, m = 500, d \\uparrow$)"))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

# Stack the two plots by row
plt.time = ggpubr::ggarrange(plt.A.time, plt.B.time, nrow = 2, ncol = 1,
                             common.legend = TRUE, legend = "none")

# Summary of the GoF measures in the two settings
plt.gof = df.gof %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  scale_colour_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  facet_grid(rows = vars(variable), cols = vars(setting),
             scale = "free", labeller = label_parsed) +
  labs(x = paste("Matrix dimension", "Latent space rank",
                 sep = paste0(rep(" ", 45), collapse = "")),
       y = "Goodness of fit measure", color = "Method") +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        strip.background = element_blank(),
        strip.text = element_text(size = 13))

# Put all together
plt = ggpubr::ggarrange(plt.time, plt.gof, nrow = 1, ncol = 2,
                        common.legend = TRUE, legend = "bottom",
                        widths = c(2,3), align = "hv", labels = "AUTO")

# Show the plot
if (SHOW) print(plt)

# Save the result
if (SAVE) {
  filepath = IMGPATH
  filename = "simulation_summary_3x3.pdf"
  zoom = 4
  nrows = 3.25
  ncols = 3
  width = 1 * ncols
  height = .6 * nrows
  ggplot2::ggsave(
    file = filename, path = filepath, plot = plt,
    width = zoom * width, height = zoom * height)
}


## END OF FILE ----

