
## Import the libraries ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(lemon)

# BiocManager::install("ComplexHeatmap")

## Load the summary data ----
{
  path <- "sim/arigoni"

  df_eigenvalues <- readRDS(file = paste0(path, "/df_eigenvalues.RDS"))
  df_model_selection <- readRDS(file = paste0(path, "/df_model_selection.RDS"))
  df_purity <- readRDS(file = paste0(path, "/df_purity.RDS"))
  df_tsne_plot <- readRDS(file = paste0(path,"/df_tsne_plot.RDS"))
  df_clust_tsne <- readRDS(file = paste0(path,"/df_clust_tsne.RDS"))
  df_tile_confusion <- readRDS(file = paste0(path,"/df_tile_confusion.RDS"))
  celltypes <- readRDS(file = paste0(path, "/celltypes.RDS"))
  plot_time <- readRDS(file = paste0(path, "/plot_time.RDS"))
  plot_heatmap <- readRDS(file = paste0(path, "/heatmap_clusters_10LF.RDS"))

}


## Plot: selection criteria ----

# Godness of fit
df_model_selection_summary <- df_model_selection %>%
  group_by(ncomp, vars) %>%
  summarise(median = median(values, na.rm = TRUE),
            ymin = min(values, na.rm = TRUE),
            ymax = max(values, na.rm = TRUE))

df_selected_points = df_model_selection_summary %>%
  filter((ncomp == 15 & vars == "AIC") |
           (ncomp == 6 & vars == "BIC") |
           (ncomp == 15 & vars == "Deviance"))

scaleFUN <- function(x) sprintf("%.1f", x)

plt_dev <- df_model_selection_summary %>%
  ggplot(map = aes(x = ncomp, y = median, ymin = ymin, ymax = ymax)) +
  geom_errorbar(width = 0.2, col = "gray40") + geom_point() +
  geom_vline(data = df_selected_points,
             map = aes(xintercept = ncomp), lty = 2, col = "red") +
  geom_errorbar(width = 0.2, col = "gray40") + geom_point() +
  geom_point(data = df_selected_points,
             map = aes(x = ncomp, y = median, col = "red")) +
  facet_grid(rows = vars(vars), scale = "free") +
  theme_bw() + xlab("Matrix rank") +
  scale_y_continuous(labels = scaleFUN) +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0.3, 0.2, 0.2, 0.3, "cm"),
        panel.spacing = unit(1, "lines"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 11),
        strip.background = element_blank())

print(plt_dev)

# Spectrum
labels_eigs <- c("Eigenvalues")
names(labels_eigs) <- "1"

plt_eigs <- df_eigenvalues %>%
  filter(ncomp %in% (c(1:10, 15, 20, 30, 40, 50))) %>%
  ggplot(map = aes(x = as.factor(ncomp), y = eigenvalues)) +
  geom_point() +
  geom_point(data = df_eigenvalues %>% filter(ncomp == 6),
             aes(x = ncomp, y = eigenvalues), col = "red") +
  geom_vline(data = df_eigenvalues %>% filter(ncomp == 6),
             aes(xintercept = ncomp), lty = 2, col = "red") +
  # geom_vline(xintercept = 6.5, col = "red") +
  facet_wrap(vars(facet), scale = "free_y",
             labeller = labeller(facet = labels_eigs),
             strip.position = "right") +
  theme_bw() +
  xlab("Matrix rank") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.2, 0.2, 0.3, "cm"),
        strip.text = element_text(size = 11),
        strip.background = element_blank())

print(plt_eigs)

# Cell-specific purity
labels_purity <- c("Mean cluster purity")
names(labels_purity) <- "1"

plt_purity <- df_purity %>%
  ggplot(map = aes(x = dim, y = mean, col = maximum)) +
  geom_point() +
  geom_line(map = aes(x = as.numeric(dim))) +
  geom_vline(xintercept = c("6", "15"), lty = 2, col = "gray40") +
  geom_point() +
  facet_wrap(~facet, scale = "free_y",
             labeller = labeller(facet = labels_purity),
             strip.position = "right") +
  labs(x = "Latent space rank", col = "Celltype") +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.2, 0.2, 0.3, "cm"),
        strip.text = element_text(size = 11),
        strip.background = element_blank())

print(plt_purity)

# Summary
plt_gof = ggpubr::ggarrange(plt_eigs, plt_dev, plt_purity,
                            nrow = 3, ncol = 1, align = "v",
                            heights = c(1,3,2))

plt_gof = ggpubr::annotate_figure(plt_gof, top = text_grob("Model selection measures", size = 14))

print(plt_gof)

## Plot: tSNE embedding ----
dim_labs_tsne <- c("Rank = 6", "Rank = 15", "Rank = 30")
names(dim_labs_tsne) <- c( "6", "15", "30")

plt_tsne_true <- df_tsne_plot %>%
  filter(dim %in% c(6, 15, 30)) %>%
  ggplot(map = aes(x  = V1, y = V2, col = celltypes)) +
  geom_point(alpha = 0.2, size = 0.5) +
  facet_grid(cols = vars(dim), labeller = labeller(dim = dim_labs_tsne)) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw() + labs(col = "Celltype") +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(plt_tsne_true)


plt_tsne_clust <- df_clust_tsne %>%
  ggplot(map = aes(x = T1, y = T2, col = membership)) +
  geom_point(alpha = 0.2, size = 0.5) +
  facet_grid(cols = vars(dims), labeller = labeller(dims = dim_labs_tsne)) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + labs(col = "Cluster") +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

print(plt_tsne_clust)

lgnd_tsne_true = g_legend(plt_tsne_true)
lgnd_tsne_clust = g_legend(plt_tsne_clust)

plt_tsne <- ggpubr::ggarrange(
  plt_tsne_true + theme(legend.position = "hiden"),
  plt_tsne_clust + theme(legend.position = "hiden"),
  nrow = 2, ncol = 1, align = "hv")

print(plt_tsne)


## Plot: confusion matrix ----
df_tile_6 = df_tile_confusion |> filter(dims == 6)
df_tile_15 = df_tile_confusion |> filter(dims == 15)
df_tile_30 = df_tile_confusion |> filter(dims == 30)


annotation_row = data.frame(Celltype = factor(unique(celltypes)))
rownames(annotation_row) = unique(celltypes)

annotation_col = data.frame(Group = factor(1:7))
rownames(annotation_col) = 1:7

mat6 = matrix(df_tile_6$hits, nrow = 7, ncol = 7)
mat15 = matrix(df_tile_15$hits, nrow = 7, ncol = 7)
mat30 = matrix(df_tile_30$hits, nrow = 7, ncol = 7)

rownames(mat6) = rownames(mat15) = rownames(mat30) = unique(celltypes)
colnames(mat6) = colnames(mat15) = colnames(mat30) = 1:7

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_dark2 = function(n) {
  RColorBrewer::brewer.pal(n = 8, name = "Dark2")
}

celltype_colors = gg_color_hue(7)
names(celltype_colors) = unique(celltypes)

group_colors = gg_color_dark2(7)
names(group_colors) = factor(1:7)

annotation_colors = list(
  Celltype = celltype_colors,
  Group = group_colors
)


myheatmap = function(mat, annotation_row, annotation_col, annotation_colors, main = "") {
  greys = grey.colors(50, start = 0.3, end = 0.95, rev = TRUE)
  mat = mat / sum(abs(mat))
  plt = pheatmap::pheatmap(mat = mat,
                           color = greys,
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           annotation_row = annotation_row,
                           annotation_col = annotation_col,
                           annotation_names_row = FALSE,
                           annotation_names_col = FALSE,
                           annotation_colors = annotation_colors,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           display_numbers = floor(100*mat),
                           number_format = "%.0f",
                           number_color = "steelblue",
                           legend = FALSE,
                           border_color = "white",
                           annotation_legend = FALSE,
                           main = main,
                           fontfamily = "serif")$gtable
  plt$grobs[[1]]$gp = grid::gpar(fontsize = 13,
                                 fontfamily = "sans",
                                 fontface = "plain")
  return (plt)
}

plt_mat6 = myheatmap(mat6, annotation_row, annotation_col, annotation_colors, main = "Rank = 6")
plt_mat15 = myheatmap(mat15, annotation_row, annotation_col, annotation_colors, main = "Rank = 15")
plt_mat30 = myheatmap(mat30, annotation_row, annotation_col, annotation_colors, main = "Rank = 30")


plt_mat = ggpubr::ggarrange(
  plt_mat6, plt_mat15, plt_mat30,
  nrow = 1, ncol = 3, common.legend = FALSE)

print(plt_mat)


## Plot: grid plot ----
plt_grid = ggpubr::ggarrange(
  ggpubr::ggarrange(plt_tsne, NULL, plt_mat, NULL,
                    nrow = 4, align = "hv", heights = c(2,0.05,1, 0.1)),
  ggpubr::ggarrange(NULL, lgnd_tsne_true, lgnd_tsne_clust, NULL,
                    nrow = 4, align = "hv", heights = c(.5,1,1,.5)),
  nrow = 1, ncol = 2, widths = c(10,2)
)

print(plt_grid)

## Plot: final layout ----
plt = ggpubr::ggarrange(plt_gof, plt_grid, ncol = 2,
                        widths = c(2,3.5), align = "h")

print(plt)

filepath = "img"
filename = "arigoni2.pdf"
zoom = 2
width = 5.5
height = 4
ggsave(filename = filename, path = filepath, plot = plt,
       width = zoom  * width, height = zoom * height)




