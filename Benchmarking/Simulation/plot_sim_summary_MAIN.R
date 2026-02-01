
## WORKSPACE SETUP ----

# Clean the workspace
rm(list = ls())
graphics.off()
options(warn = -1)

# Packages
library(ggplot2)
library(ggrepel)
library(dplyr)
library(latex2exp)

# Default theme
theme_set(theme_bw())

# Global variables
SAVE = TRUE
SHOW = TRUE

# File and image paths
FILEPATH = paste("Benchmarking", "Simulation", "data", sep = "/")
IMGPATH = paste("Benchmarking", "Simulation", "img", sep = "/")

# Color palette
COLORS = hcl(h = seq(15, 375, length = 13), l = 65, c = 100)[1:12]

## UTILITY FUNCTIONS ----

# Loads the data corresponding to:
# - n: data matrix rows
# - m: data matrix columns
# - d: fixed latent space rank
# - i: number of iterations in the simulation experiment
load.data = function (n, m, d, i) {

  # Define the file path and name
  filepath = FILEPATH
  filezip = "summary_sim_MAIN.zip"
  fileid = paste("_n", n, "_m", m, "_d", d, "_i", i, sep = "")
  filename = paste("summary_sim", fileid, ".csv", sep = "")

  # Load the data-set
  df = read.table(unz(paste(filepath, filezip, sep = "/"), filename),
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
  old.names = c("Model", "Time", "Memory", "RSS", "Dev", "Sil", "Purity", "Iter", "Dim", "NComp")
  new.names = c("model", "time", "memory", "error", "deviance", "silhouette", "purity", "iteration", "dimension", "rank")
  old.models = c("CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "GFM-AM", "GFM-VEM", "COAP", "AIRWLS", "Newton", "SGD")
  new.models = c("CMF", "NMF", "NMF+", "AvaGrad", "Fisher", "NBWaVE", "GFM-AM", "GFM-VEM", "COAP", "AIRWLS", "Newton", "aSGD")

  # Filter out the unused rows and columns
  # (do not change the order of the following lines of code)
  # df = df[-which(df$Set != "Test"),]
  # df = df[-which(df$Model == "C-SGD"),]
  # df = df[,-which(colnames(df) == "Set")]
  # df = df[,-which(colnames(df) == "Cos")]
  df = df %>% dplyr::filter(Set == "Test")
  df = df %>% dplyr::filter(Model != "C-SGD")
  df = df %>% dplyr::select(-Set, -Cos)
  df = df %>% dplyr::mutate(Memory = Memory * 1.048576)

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

  # Fill extreme values with NAs
  df$error[df$error>1.5] = NA
  df$deviance[df$deviance>1.5] = NA


  # Rebuild the data-frame changing the order of the columns
  # df = data.frame(
  #   "iteration" = df$iteration,
  #   "model" = df$model,
  #   "dimension" = df$dimension,
  #   "rank" = df$rank,
  #   "time" = df$time,
  #   "memory" = df$memory,
  #   "error" = df$error,
  #   "deviance" = df$deviance,
  #   "silhouette" = df$silhouette,
  #   "purity" = df$purity)
  df = df %>% dplyr::select(iteration, model, dimension, rank, time,
                            memory, error, deviance, silhouette, purity)

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
              memory = mean(memory, na.rm = TRUE),
              error = mean(100*error, na.rm = TRUE),
              deviance = mean(100*deviance, na.rm = TRUE),
              silhouette = mean(silhouette, na.rm = TRUE),
              purity = mean(purity, na.rm = TRUE),
              .groups = "drop") %>%
    as.data.frame()

  # Do you want to melt the dataset
  if (melt) {
    idxA = which(df$setting == "A")
    idxB = which(df$setting == "B")
    vars = c("Time", "Memory", "Error", "Deviance", "Silhouette", "Purity")
    dims = rep(NA, times = nrow(df))
    dims[idxA] = df$dimension[idxA]
    dims[idxB] = df$rank[idxB]

    df$dimension = dims
    df = df %>% dplyr::select(-rank)
    # df = df %>% select(-error)

    df = reshape2::melt(df, id.vars = c("model", "setting", "dimension"))
  }

  return (df)
}

## DATA LOAD ----

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


## PLOT: DATA EXPLORATION ----

# Merge the two simulation settings
df.flat.tmp = df.flat %>% mutate(error = 5*error*(dimension)^-.4)
df.deep.tmp = df.deep %>% mutate(error = 5*error*(dimension)^-.4)
df.melt = merge.data(df.flat.tmp, df.deep.tmp, melt = TRUE)


# Map the variable labels to upper case letters
df.melt$variable = factor(df.melt$variable,
                          levels = c("time", "memory", "error", "deviance", "silhouette", "purity"),
                          labels = c("Time", "Memory", "Error", "Deviance", "Silhouette", "Purity"))

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


## DATA PREPARATION ----

# Separate the time summary to the other gof measures
df.time = df.melt %>%
  filter(variable == "Time") %>%
  mutate(label = if_else(setting == levels(setting)[1] & dimension != 100, NA, model)) %>%
  mutate(label = if_else(setting == levels(setting)[2] & dimension != 25, NA, model))

df.memory = df.melt %>%
  filter(variable == "Memory") %>%
  mutate(label = if_else(setting == levels(setting)[1] & dimension != 100, NA, model)) %>%
  mutate(label = if_else(setting == levels(setting)[2] & dimension != 25, NA, model))

df.gof = df.melt %>%
  filter(variable != "Time", variable != "Memory") %>%
  mutate(label = if_else(setting == levels(setting)[1] & dimension != 100, NA, model)) %>%
  mutate(label = if_else(setting == levels(setting)[2] & dimension != 25, NA, model))

# Separate the execution times of the two settings
df.A.time = df.melt %>% filter(variable == "Time", setting == levels(setting)[1])
df.B.time = df.melt %>% filter(variable == "Time", setting == levels(setting)[2])

df.A.time$label = as.character(df.A.time$model)
df.B.time$label = as.character(df.B.time$model)
df.A.time$label[df.A.time$dimension != 100] <- NA
df.B.time$label[df.B.time$dimension != 25] <- NA

# Separate the memory consumption of the two settings
df.A.memory = df.melt %>% filter(variable == "Memory", setting == levels(setting)[1])
df.B.memory = df.melt %>% filter(variable == "Memory", setting == levels(setting)[2])

df.A.memory$label = as.character(df.A.memory$model)
df.B.memory$label = as.character(df.B.memory$model)
df.A.memory$label[df.A.memory$dimension != 100] <- NA
df.B.memory$label[df.B.memory$dimension != 25] <- NA

# Separate the gof of the two settings
df.A.gof = df.melt %>% filter(variable != "Time", variable != "Memory", setting == levels(setting)[1])
df.B.gof = df.melt %>% filter(variable != "Time", variable != "Memory", setting == levels(setting)[2])

df.A.gof$label = as.character(df.A.gof$model)
df.B.gof$label = as.character(df.B.gof$model)
df.A.gof$label[df.A.gof$dimension != 100] <- NA
df.B.gof$label[df.B.gof$dimension != 25] <- NA

## PLOT: TIME & MEMORY ----

df.A.time_memory = rbind(df.A.time, df.A.memory) %>%
  mutate(variable = factor(variable, levels = c("Time", "Memory"),
                           labels = c("Execution time (s)",
                                      "Peak memory consumption (MB)")))

df.B.time_memory = rbind(df.B.time, df.B.memory) %>%
  mutate(variable = factor(variable, levels = c("Time", "Memory"),
                           labels = c("Execution time (s)",
                                      "Peak memory consumption (MB)")))

# Execution time in the first setting
plt.A.time_memory = df.A.time_memory %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  geom_text_repel(map = aes(label = label), hjust=-.2,
                  direction="y", show.legend = FALSE) +
  facet_grid(rows = vars(variable), scales = "free") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25)),
                     breaks = c(10, 25, 50, 75, 100)) +
  scale_colour_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(x = "Matrix dimension", y = "Execution time (s)", color = "Method",
       title = TeX(paste("Setting A ($n \\uparrow, m \\uparrow, d = 5$)"))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = "transparent"),
        strip.text = element_text(size = 13))

plt.B.time_memory = df.B.time_memory %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  geom_text_repel(map = aes(label = label), hjust=-.2,
                  direction="y", show.legend = FALSE) +
  facet_grid(rows = vars(variable), scales = "free") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25)),
                     breaks = c(5, 10, 15, 20, 25)) +
  scale_colour_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) +
  labs(x = "Latent space rank", y = "Execution time (s)", color = "Method",
       title = TeX(paste("Setting B ($n = 5000, m = 500, d \\uparrow$)"))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = "transparent"),
        strip.text = element_text(size = 13))

plt = ggpubr::ggarrange(plt.A.time_memory, plt.B.time_memory, nrow = 1, ncol = 2,
                        common.legend = TRUE, legend = "none",
                        align = "hv", labels = "AUTO")


# Show the plot
if (SHOW) print(plt)

# Save the result
if (SAVE) {
  filepath = IMGPATH
  filename = "simulation_time_memory.pdf"
  zoom = 3.75
  nrows = 3
  ncols = 3
  width = 1 * ncols
  height = .6 * nrows
  ggplot2::ggsave(
    file = filename, path = filepath, plot = plt,
    width = zoom * width, height = zoom * height)
}


## PLOT: GOODNESS OF FIT ----

# Summary of the GoF measures in the two settings
plt.A.gof = df.A.gof %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  facet_grid(rows = vars(variable),
             scale = "free", labeller = label_parsed) +
  geom_text_repel(map = aes(label = label), hjust=-.2,
                  direction="y", show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25)),
                     breaks = c(10, 25, 50, 75, 100)) +
  labs(x = "Matrix dimension", y = "Goodness of fit measure", color = "Method",
       title = TeX(paste("Setting A ($n \\uparrow, m \\uparrow, d = 5$)"))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = "transparent"),
        strip.text = element_text(size = 13))

plt.B.gof = df.B.gof %>%
  ggplot(map = aes(x = dimension, y = value, color = model)) +
  geom_line() + geom_point(size = 3) +
  facet_grid(rows = vars(variable),
             scale = "free", labeller = label_parsed) +
  geom_text_repel(map = aes(label = label), hjust=-.2,
                  direction="y", show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25)),
                     breaks = c(5, 10, 15, 20, 25)) +
  labs(x = "Latent space rank", y = "Goodness of fit measure", color = "Method",
       title = TeX(paste("Setting B ($n = 5000, m = 500, d \\uparrow$)"))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 13),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = "transparent"),
        strip.text = element_text(size = 13))

plt = ggpubr::ggarrange(plt.A.gof, plt.B.gof,
                        nrow = 1, ncol = 2, common.legend = TRUE,
                        legend = "none", align = "hv", labels = "AUTO")


# Show the plot
if (SHOW) print(plt)

# Save the result
if (SAVE) {
  filepath = IMGPATH
  filename = "simulation_gof_summary.pdf"
  zoom = 3.75
  nrows = 5
  ncols = 3
  width = 1 * ncols
  height = .6 * nrows
  ggplot2::ggsave(
    file = filename, path = filepath, plot = plt,
    width = zoom * width, height = zoom * height)
}



## END OF FILE ----
