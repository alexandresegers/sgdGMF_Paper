#' ---
#' title: "Untitled"
#' output: html_document
#' date: '2024-04-17'
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## -------------------------------------------------------------------------------------------------
deviance <- list()
time <- list()
for(i in c("glmpca_100k_avagrad_def_1", "glmpca_100k_avagrad_def_2", 
           "glmpca_100k_avagrad_def_3", "glmpca_100k_avagrad_def_4", 
           "glmpca_100k_avagrad_def_5", 
           "glmpca_200k_avagrad_def_1", "glmpca_200k_avagrad_def_2", 
           "glmpca_200k_avagrad_def_3", "glmpca_200k_avagrad_def_4", 
           "glmpca_200k_avagrad_def_5", 
           "glmpca_300k_avagrad_def_1", "glmpca_300k_avagrad_def_2", 
           "glmpca_300k_avagrad_def_3", "glmpca_300k_avagrad_def_4", 
           "glmpca_300k_avagrad_def_5", 
           "glmpca_100k_fisher_def_1", "glmpca_100k_fisher_def_2", 
           "glmpca_100k_fisher_def_3", "glmpca_100k_fisher_def_4", 
           "glmpca_100k_fisher_def_5", 
           "glmpca_200k_fisher_def_1", "glmpca_200k_fisher_def_2", 
           "glmpca_200k_fisher_def_3", "glmpca_200k_fisher_def_4", 
           "glmpca_200k_fisher_def_5", 
           "glmpca_300k_fisher_def_1", "glmpca_300k_fisher_def_2", 
           "glmpca_300k_fisher_def_3", "glmpca_300k_fisher_def_4", 
           "glmpca_300k_fisher_def_5",
           "sgdgmffit_100k_samples_1000_maxiter_1000_250_size_1",
           "sgdgmffit_100k_samples_1000_maxiter_1000_250_size_2",
           "sgdgmffit_100k_samples_1000_maxiter_1000_250_size_3",
           "sgdgmffit_100k_samples_1000_maxiter_1000_250_size_4",
           "sgdgmffit_100k_samples_1000_maxiter_1000_250_size_5",
           "sgdgmffit_200k_samples_2000_maxiter_1000_250_size_1",
           "sgdgmffit_200k_samples_2000_maxiter_1000_250_size_2",
           "sgdgmffit_200k_samples_2000_maxiter_1000_250_size_3",
           "sgdgmffit_200k_samples_2000_maxiter_1000_250_size_4",
           "sgdgmffit_200k_samples_2000_maxiter_1000_250_size_5",
           "sgdgmffit_300k_samples_3000_maxiter_1000_250_size_1",
           "sgdgmffit_300k_samples_3000_maxiter_1000_250_size_2",
           "sgdgmffit_300k_samples_3000_maxiter_1000_250_size_3",
           "sgdgmffit_300k_samples_3000_maxiter_1000_250_size_4",
           "sgdgmffit_300k_samples_3000_maxiter_1000_250_size_5")){
  fit <- readRDS(file = paste0("Output_files/",i, ".RDS"))
  deviance[[paste0(i)]] <- fit$fit$dev
  time[[paste0(i)]] <- fit$time
}
for(i in c(
           "NewWave_100k_def", "NewWave_100k_def_2",
           "NewWave_100k_def_3", "NewWave_100k_def_4",
           "NewWave_100k_def_5",
           "NewWave_200k_def", "NewWave_200k_def_2",
           "NewWave_200k_def_3", "NewWave_200k_def_4",
           "NewWave_200k_def_5",
           "NewWave_300k_def", "NewWave_300k_def_2",
           "NewWave_300k_def_3", "NewWave_300k_def_4",
           "NewWave_300k_def_5")){
  fit <- readRDS(file = paste0("Output_files/", i, ".RDS"))
  time[[paste0(i)]] <- fit$time
}
saveRDS(deviance, file = "Output_files/deviances_comparison_samples.RDS")
saveRDS(time, file = "Output_files/time_comparison_samples.RDS")


#' 
## -------------------------------------------------------------------------------------------------
df_time <- data.frame("time" = c(time$glmpca_100k_avagrad_def_1[3]/60, 
                                 time$glmpca_100k_avagrad_def_2[3]/60, 
                                 time$glmpca_100k_avagrad_def_3[3]/60, 
                                 time$glmpca_100k_avagrad_def_4[3]/60, 
                                 time$glmpca_100k_avagrad_def_5[3]/60, 
                                 time$glmpca_200k_avagrad_def_1[3]/60, 
                                 time$glmpca_200k_avagrad_def_2[3]/60,
                                 time$glmpca_200k_avagrad_def_3[3]/60,
                                 time$glmpca_200k_avagrad_def_4[3]/60,
                                 time$glmpca_200k_avagrad_def_5[3]/60,
                                 time$glmpca_300k_avagrad_def_1[3]/60,
                                 time$glmpca_300k_avagrad_def_2[3]/60,
                                 time$glmpca_300k_avagrad_def_3[3]/60,
                                 time$glmpca_300k_avagrad_def_4[3]/60,
                                 time$glmpca_300k_avagrad_def_5[3]/60,
                                 time$glmpca_100k_fisher_def_1[3]/60,
                                 time$glmpca_100k_fisher_def_2[3]/60,
                                 time$glmpca_100k_fisher_def_3[3]/60,
                                 time$glmpca_100k_fisher_def_4[3]/60,
                                 time$glmpca_100k_fisher_def_5[3]/60,
                                 time$glmpca_200k_fisher_def_1[3]/60,
                                 time$glmpca_200k_fisher_def_2[3]/60,
                                 time$glmpca_200k_fisher_def_3[3]/60,
                                 time$glmpca_200k_fisher_def_4[3]/60,
                                 time$glmpca_200k_fisher_def_5[3]/60,
                                 time$glmpca_300k_fisher_def_1[3]/60,
                                 time$glmpca_300k_fisher_def_2[3]/60,
                                 time$glmpca_300k_fisher_def_3[3]/60,
                                 time$glmpca_300k_fisher_def_4[3]/60,
                                 time$glmpca_300k_fisher_def_5[3]/60,
                                 time$NewWave_100k_def[3]/60,
                                 time$NewWave_100k_def_2[3]/60,
                                 time$NewWave_100k_def_3[3]/60,
                                 time$NewWave_100k_def_4[3]/60,
                                 time$NewWave_100k_def_5[3]/60,
                                 time$NewWave_200k_def[3]/60,
                                 time$NewWave_200k_def_2[3]/60,
                                 time$NewWave_200k_def_3[3]/60,
                                 time$NewWave_200k_def_4[3]/60,
                                 time$NewWave_200k_def_5[3]/60,
                                 time$NewWave_300k_def[3]/60,
                                 time$NewWave_300k_def_2[3]/60,
                                 time$NewWave_300k_def_3[3]/60,
                                 time$NewWave_300k_def_4[3]/60,
                                 time$NewWave_300k_def_5[3]/60,
                                 time$sgdgmffit_100k_samples_1000_maxiter_1000_250_size_1,
                                 time$sgdgmffit_100k_samples_1000_maxiter_1000_250_size_2,
                                 time$sgdgmffit_100k_samples_1000_maxiter_1000_250_size_3,
                                 time$sgdgmffit_100k_samples_1000_maxiter_1000_250_size_4,
                                 time$sgdgmffit_100k_samples_1000_maxiter_1000_250_size_5,
                                 time$sgdgmffit_200k_samples_2000_maxiter_1000_250_size_1,
                                 time$sgdgmffit_200k_samples_2000_maxiter_1000_250_size_2,
                                 time$sgdgmffit_200k_samples_2000_maxiter_1000_250_size_3,
                                 time$sgdgmffit_200k_samples_2000_maxiter_1000_250_size_4,
                                 time$sgdgmffit_200k_samples_2000_maxiter_1000_250_size_5,
                                 time$sgdgmffit_300k_samples_3000_maxiter_1000_250_size_1,
                                 time$sgdgmffit_300k_samples_3000_maxiter_1000_250_size_2,
                                 time$sgdgmffit_300k_samples_3000_maxiter_1000_250_size_3,
                                 time$sgdgmffit_300k_samples_3000_maxiter_1000_250_size_4,
                                 time$sgdgmffit_300k_samples_3000_maxiter_1000_250_size_5),
                      "samples" = rep(rep(c(100000, 200000, 300000), each = 5), times = 4),
                      "method" = factor(rep(c("Avagrad", "Fisher", "NewWave", "SGD"), each = 5*3)))

#' 
#' 
## -------------------------------------------------------------------------------------------------
library(dplyr)
median_time <- df_time %>% group_by(samples, method) %>% summarise("median" = median(time))
library(ggplot2)
plot_time_1 <- ggplot(df_time, aes(x = samples, y = time, col = method)) + 
  geom_point() + 
  geom_line(data = median_time, aes(x = as.numeric(samples), y = median, col = method)) + theme_bw() + 
  xlab("Number of cells") + 
  ylab("Time (min)") + 
  labs(col = "Analysis method")
ggsave(plot_time_1, filename = "Figures/plot_time_1.pdf")
saveRDS(plot_time_1, file = "Figures/plot_time_1.RDS")
  
plot_time_2 <- ggplot(df_time, aes(x = samples, y = log10(time), col = method)) + 
  geom_point() + 
  geom_line(data = median_time, aes(x = as.numeric(samples), y = log10(median), col = method)) + theme_bw() + 
  xlab("Number of cells") + 
  ylab("log10(Time) (min)") + 
  labs(col = "Analysis method")
ggsave(plot_time_2, filename = "Figures/plot_time_2.pdf")
saveRDS(plot_time_2, file = "Figures/plot_time_2.RDS")


#' 
