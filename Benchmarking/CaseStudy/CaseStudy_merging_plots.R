#' ---
#' title: "Untitled"
#' output: html_document
#' date: '2024-04-25'
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## -------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

#' 
## -------------------------------------------------------------------------------------------------
eigenvalues <- readRDS(file = "Output_files/eigenvalues_sgdGMF_link.RDS")


## -------------------------------------------------------------------------------------------------
p1 <- ggplot(data = data.frame(eigenvalues = eigenvalues$lambda[1:50]), 
             aes(x = 1:50, y = eigenvalues)) + 
    geom_point() + 
    theme_bw() + 
    ylab("Eigenvalue") + 
    xlab("Component number") + 
    theme(panel.grid = element_blank()) + 
    #geom_vline(xintercept = 10.5, col = "red") + 
    geom_segment(
               aes(x = 10.5, 
                   xend  = 10.5, 
                   y = 0, 
                   yend = 30),
               color = "red", lwd = 0.5) + 

    geom_rect(aes(xmin = 2.5, xmax = 50.5, ymin = 0, ymax = 30), color = "grey", alpha = 0,
              linetype='dashed', linewidth = 0.2) 
p1

p2 <- ggplot(data = data.frame(eigenvalues = eigenvalues$lambda[3:50]),
             aes(x = 3:50, y = eigenvalues)) + geom_point() +
    theme_bw() + 
    ylab("Eigenvalue") + 
    xlab("Component number") + 
    theme(panel.grid = element_blank()) + 
    geom_vline(xintercept = 10.5, col = "red")

p2
p3 <- p1  + annotation_custom(ggplotGrob(p2), xmin = 20, xmax = 50.5, ymin = 70, ymax = 110) +
  geom_rect(aes(xmin = 20, xmax = 50.5, ymin = 70, ymax = 110), color='grey', 
            linetype='dashed', alpha=0, linewidth = 0.2) +
  geom_path(aes(x,y,group=grp), 
            data=data.frame(x = c(2.5,20,50.5,50.5), y=c(30,70,27,70),grp=c(1,1,2,2)),
            linetype='dashed', color = "grey", linewidth = 0.2)

ggsave(p3, filename = "Figures/scree_plot.pdf")


#' 
## -------------------------------------------------------------------------------------------------

plot_cluster <- readRDS(file = "Figures/heatmap_clusters_10LF_2.RDS")
plot_time <- readRDS(file = "Figures/plot_time_1.RDS")

plot_casestudy <- ggarrange(ggplotify::as.ggplot(plot_cluster), 
                            plot_time + theme(legend.position = "top",
                                              legend.text = element_text(size = 10),
                                              legend.title = element_text(size = 12),
                                              panel.grid = element_blank(),
                                              plot.margin = margin(r = 10)) + 
    guides(color=guide_legend(nrow=2, byrow=TRUE, title.position ="top", title.hjust = 0.5)), 
                            ncol = 2, 
                            widths = c(3,2),labels = "AUTO")

ggsave(plot_casestudy, filename = "Figures/CaseStudy_mainfigure.pdf")

#' 
