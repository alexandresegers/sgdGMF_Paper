# List all files in the current directory

dir <- "~/sgdGMF_Paper/Benchmarking/CaseStudy/Memory/Output_files/"
all_files <- list.files(path = dir, full.names = FALSE)
size <- c(100, 200, 300)

results_avagrad <- matrix(NA, ncol = length(size), nrow= 5)
time_avagrad <- matrix(NA, ncol = length(size), nrow= 5)

for(i in seq_len(length(size))){
  pattern <- paste0("^glmpca_", size[i], "k_avagrad_def_[0-9]+\\.RDS$")
  matching_files <- grep(pattern, all_files, value = TRUE)

  for(j in seq_len(length(matching_files))){
    lines <- readRDS(paste0(dir, matching_files[j]))
    rss_line <- lines$memory[[4]]
    time <- lines$time[3]

    results_avagrad[j,i] <- rss_line
    time_avagrad[j,i] <- time
  }
}


results_fisher <- matrix(NA, ncol = length(size), nrow= 5)
time_fisher <- matrix(NA, ncol = length(size), nrow= 5)

for(i in seq_len(length(size))){
  pattern <- paste0("^glmpca_", size[i], "k_fisher_def_[0-9]+\\.RDS$")
  matching_files <- grep(pattern, all_files, value = TRUE)

  for(j in seq_len(length(matching_files))){
    lines <- readRDS(paste0(dir, matching_files[j]))
    rss_line <- lines$memory[[4]]
    time <- lines$time[3]

    results_fisher[j,i] <- rss_line
    time_fisher[j,i] <- time
  }
}

results_nbwave <- matrix(NA, ncol = length(size), nrow= 5)
time_nbwave <- matrix(NA, ncol = length(size), nrow= 5)

for(i in seq_len(length(size))){
  pattern <- paste0("^NewWave_", size[i], "k_def_[0-9]+\\.RDS$")
  matching_files <- grep(pattern, all_files, value = TRUE)

  for(j in seq_len(length(matching_files))){
    lines <- readRDS(paste0(dir, matching_files[j]))
    rss_line <- lines$memory[[4]]
    time <- lines$time[3]

    results_nbwave[j,i] <- rss_line
    time_nbwave[j,i] <- time

  }
}


results_sgdGMF <- matrix(NA, ncol = length(size), nrow= 5)
time_sgdGMF <- matrix(NA, ncol = length(size), nrow= 5)

for(i in seq_len(length(size))){
  pattern <- paste0("^sgdgmffit_", size[i], "k_memory[0-9]+\\.RDS$")
  matching_files <- grep(pattern, all_files, value = TRUE)

  for(j in seq_len(length(matching_files))){
    lines <- readRDS(paste0(dir, matching_files[j]))
    rss_line <- lines$memory[[4]]
    time <- lines$time

    results_sgdGMF[j,i] <- rss_line
    time_sgdGMF[j,i] <- time
  }
}




memory_usage <- matrix(as.numeric(sapply(FUN = rbind, X = list(results_sgdGMF, results_avagrad, results_fisher, results_nbwave)))*1.048576/1024, ncol = 4)
colnames(memory_usage) <- c("sgdGMF", "Avagrad", "Fisher", "NBWaVe")
memory_usage


time_usage <- matrix(as.numeric(sapply(FUN = rbind, X = list(time_sgdGMF*60,
                                                               time_avagrad,
                                                               time_fisher,
                                                               time_nbwave)))/60,
                       ncol = 4)
colnames(time_usage) <- c("sgdGMF", "Avagrad", "Fisher", "NBWaVe")
time_usage


df_memory_usage <- data.frame("Memory" = c(memory_usage),
                              "Method" = rep(c("aSGD", "AvaGrad", "Fisher", "NBWaVe"), each = 5*3),
                              "Samples" = (rep(rep(c(100000, 200000, 300000), each = 5), times = 4))
)

df_time_usage <- data.frame("Time" = c(time_usage),
                              "Method" = rep(c("aSGD", "AvaGrad", "Fisher", "NBWaVe"), each = 5*3),
                              "Samples" = (rep(rep(c(100000, 200000, 300000), each = 5), times = 4))
)


df_memory_usage$Method <- factor(df_memory_usage$Method, levels = c( "AvaGrad", "Fisher", "NBWaVe", "aSGD"))
df_time_usage$Method <- factor(df_time_usage$Method, levels = c( "AvaGrad", "Fisher", "NBWaVe", "aSGD"))

library(RColorBrewer)
library(dplyr)

colScale <- scale_colour_manual(name = "Method",
                                values = c("AvaGrad" = "#7CAE00",
                                           "Fisher" = "#00BA38",
                                           "NBWaVe" = "#00C08B",
                                           "COAP" = "#619CFF",
                                           "aSGD" = "#FF64B0"))

library(ggplot2)
library(dplyr)

median_memory <- df_memory_usage %>% group_by(Samples, Method) %>% summarise("median" = median(Memory))
median_time <- df_time_usage %>% group_by(Samples, Method) %>% summarise("median" = median(Time))


plot <- ggplot(data = df_memory_usage, aes(x = Samples, y = Memory, col = Method, group = interaction(Method))) +
  geom_point() +
  geom_line(data = median_memory, aes(x = as.numeric(Samples), y = median, col = Method)) + theme_bw() +
  colScale +
  xlab("Number of cells") +
  ylab("Memory (GB)") +
  labs(col = "Analysis method")

plot_time <- ggplot(df_time_usage, aes(x = Samples, y = Time, col = Method, group = interaction(Method))) +
  geom_point() +
  geom_line(data = median_time, aes(x = as.numeric(Samples), y = median, col = Method)) + theme_bw() +
  colScale +
  xlab("Number of cells") +
  ylab("Time (min)") +
  labs(col = "Analysis method")

arranged_plot <- ggarrange(plot + theme(legend.position = "none"),
          plot_time + theme(legend.position = "none"),
          get_legend(plot + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))),
          nrow = 3,
          heights = c(0.45, 0.45, 0.1))

saveRDS(df_memory_usage, "~/sgdGMF_Paper/Benchmarking/CaseStudy/Memory/df_memory.RDS")
saveRDS(df_time_usage, "~/sgdGMF_Paper/Benchmarking/CaseStudy/Memory/df_time.RDS")

saveRDS(arranged_plot, "~/sgdGMF_Paper/Benchmarking/CaseStudy/Memory/time_memory_plot.RDS")

ggsave(arranged_plot, filename = "~/sgdGMF_Paper/Benchmarking/CaseStudy/Memory/plot_memory.pdf")

