
# plot_rsem_sim_benchmark.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/plots/mapping/")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_tsv(filename)
  
  if (dir_split[2] == "vg") {
  
    data <- data %>%
      add_column(Reads = dir_split[5]) %>%
      add_column(Method = dir_split[6]) %>%
      add_column(Graph = dir_split[7])

  } else {
    
    data <- data %>%
      add_column(Method = dir_split[2]) %>%
      add_column(Reads = dir_split[5]) %>%
      add_column(Graph = dir_split[6])
  }

  data$Graph = recode_factor(data$Graph, "gencode100" = "Linear", "1kg_NA12878_gencode100" = "Personal", "1kg_nonCEU_af001_gencode85" = "1000g (85% sj)", "1kg_nonCEU_af001_gencode100" = "1000g")
  
  return(data)
}

wes_cols <- c(wes_palette("Darjeeling2"), wes_palette("Darjeeling1"))

overlap_data_o1 <- map_dfr(list.files(pattern=".*_stats_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o1$Threshold <- "Overlap >= 1%"

overlap_data_o1 <- overlap_data_o1 %>%
  mutate(Correct = Overlap > 0.01)

overlap_data_o50 <- map_dfr(list.files(pattern=".*_stats_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o50$Threshold <- "Overlap >= 50%"

overlap_data_o50 <- overlap_data_o50 %>%
  mutate(Correct = Overlap > 0.5)

overlap_data_o90 <- map_dfr(list.files(pattern=".*_stats_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o90$Threshold <- "Overlap >= 90%"

overlap_data_o90 <- overlap_data_o90 %>%
  mutate(Correct = Overlap > 0.9)

overlap_data_o99 <- map_dfr(list.files(pattern=".*_stats_rsem_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o99$Threshold <- "Overlap >= 99%"

overlap_data_o99 <- overlap_data_o99 %>%
  mutate(Correct = Overlap >= 0.99)



overlap_data <- bind_rows(overlap_data_o1, overlap_data_o50, overlap_data_o90, overlap_data_o99)

overlap_data$Threshold <- as.factor(overlap_data$Threshold)

overlap_data[overlap_data$MapQ == 255,]$MapQ <- 60

overlap_data_methods <- overlap_data %>%
  filter(Method != "mpmap")  %>%
  filter(Method != "mpmap_old") %>%
  filter(Method != "mpmap_fast") %>%
  filter(Method != "mpmap_fast2") %>%
  filter(Method != "mpmap_fast3") %>%
  filter(Method != "mpmap_fast4") %>%
  filter(Method != "mpmap_fast6") %>%
  filter(Method != "map") 

overlap_data_methods$Method = recode_factor(overlap_data_methods$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap_fast7" = "vg mpmap", "mpmap_final1" = "vg mpmap (final1)")

plotOverlapBenchmark(overlap_data_methods, wes_cols, "rsem_sim_benchmark_overlap_methods.pdf")

