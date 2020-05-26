
# plot_vg_sim_benchmark.R

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
  
  data <- read.table(filename, header = F)
  colnames(data) <- c("Count", "Overlap", "MapQ", "Method")
  
  if (dir_split[2] == "vg") {
    
    data <- data %>%
      add_column(Reads = dir_split[4]) %>%
      add_column(Graph = dir_split[6])
    
  } else {
    
    data <- data %>%
      add_column(Reads = dir_split[4]) %>%
      add_column(Graph = dir_split[5])
  }
  
  data$Method = recode_factor(data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap", "mpmap_old" = "vg mpmap (old)", "mpmap_fast" = "vg mpmap (fast)", "mpmap_fast2" = "vg mpmap (fast2)", "mpmap_fast3" = "vg mpmap (fast3)", "mpmap_fast4" = "vg mpmap (fast4)", "mpmap_fast7" = "vg mpmap (fast7)")
  data$Graph = recode_factor(data$Graph, "gencode100" = "Linear", "gencode100_genes" = "Splice", "1kg_nonCEU_af001_gencode100" = "1000g")
  
  return(data)
}

wes_cols <- c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2"))

distance_data_d100 <- map_dfr(list.files(pattern=".*_vg_d100_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d100$Threshold <- "Distance <= 100"

distance_data_d10 <- map_dfr(list.files(pattern=".*_vg_d10_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d10$Threshold <- "Distance <= 10"

distance_data <- rbind(distance_data_d100, distance_data_d10)
distance_data$Threshold <- as.factor(distance_data$Threshold)

distance_data <- distance_data %>%
  mutate(Correct = Overlap > 0.5) 

distance_data_methods <- distance_data %>%
  filter(Method != "STAR") %>%
  filter(Graph == "1000g")

plotDistanceBenchmark(distance_data_methods, wes_cols, "vg_sim_benchmark_distance_methods.pdf")

overlap_data_graphs <- distance_data %>%
  filter(Method == "HISAT2" | Method == "vg mpmap (fast4)") %>%
  filter(Threshold == "Distance <= 10")

plotDistanceBenchmark(overlap_data_graphs, wes_cols, "vg_sim_benchmark_distance_graphs.pdf")

