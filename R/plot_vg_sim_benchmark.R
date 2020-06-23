
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
  
  data <- data %>%
    add_column(Reads = dir_split[5])  
  
  if (dir_split[2] == "vg") {
    
    data <- data %>%
      add_column(Param = "") %>%
      add_column(Graph = dir_split[7])

  } else if (dir_split[5] == "ERR187607_uni" | dir_split[5] == "ENCSR958UOC_rep1_uni") {
  
    data <- data %>%
      add_column(Param = dir_split[6]) %>%
      add_column(Graph = dir_split[7])
    
  } else {
    
    data <- data %>%
      add_column(Param = "") %>%
      add_column(Graph = dir_split[6])
  }
  
  data <- data %>%
    rowwise() %>%
    mutate(Method = paste(Method, Param, collapse = " "))
  
  #data$Method = recode_factor(data$Method, "hisat2" = "HISAT2", "star" = "STAR", "bowtie2" = "Bowtie2", "map" = "vg map", "mpmap" = "vg mpmap", "mpmap_old" = "vg mpmap (old)", "mpmap_fast" = "vg mpmap (fast)", "mpmap_fast2" = "vg mpmap (fast2)", "mpmap_fast3" = "vg mpmap (fast3)", "mpmap_fast4" = "vg mpmap (fast4)", "mpmap_fast7" = "vg mpmap (fast7)")
  data$Graph = recode_factor(data$Graph, "gencode100" = "Linear", "linear" = "Linear", "gencode100_genes" = "Splice", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (non-CEU)")
  
  return(data)
}

wes_cols <- c(wes_palette("Darjeeling2"), wes_palette("Darjeeling1"))

distance_data_d100 <- map_dfr(list.files(pattern=".*_vg_d100_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d100$Threshold <- "Distance <= 100"

distance_data_d20 <- map_dfr(list.files(pattern=".*_vg_d20_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d20$Threshold <- "Distance <= 20"

distance_data_d10 <- map_dfr(list.files(pattern=".*_vg_d10_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d10$Threshold <- "Distance <= 10"

distance_data_d1 <- map_dfr(list.files(pattern=".*_vg_d1_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d1$Threshold <- "Distance <= 1"

#distance_data <- rbind(distance_data_d10, distance_data_d1)
#distance_data <- rbind(distance_data, distance_data_d1)

distance_data <- bind_rows(distance_data_d1, distance_data_d10, distance_data_d20, distance_data_d100) %>%
  mutate(Correct = Overlap > 0.5) 

distance_data$Threshold <- as.factor(distance_data$Threshold)

distance_data_polya <- distance_data %>%
  filter(Reads == "sim_SRR1153470_uni_vg") %>%
  filter(Method != "mpmap ")  %>%
  filter(Method != "mpmap_old ") %>%
  filter(Method != "mpmap_fast ") %>%
  filter(Method != "mpmap_fast2 ") %>%
  filter(Method != "mpmap_fast3 ") %>%
  filter(Method != "mpmap_fast4 ") %>%
  filter(Method != "mpmap_fast6 ") %>%
  filter(Method != "star ")  

distance_data_polya$Method = recode_factor(distance_data_polya$Method, "hisat2 " = "HISAT2", "star " = "STAR", "map " = "vg map", "mpmap_fast7 " = "vg mpmap", "mpmap_final1 " = "vg mpmap (final1)")

plotDistanceBenchmark(distance_data_polya, wes_cols, "vg_sim_benchmark_distance.pdf")


distance_data_mirna <- distance_data %>%
  filter(Method != "bowtie2 vse") %>%
  filter(Method != "bowtie2 vsl") %>%
  filter(Method != "map map") %>%
  filter(Method != "star ns") %>%
  filter(Method != "mpmap ") %>%
  filter(Method != "mpmap_seed ")

distance_data_mirna$Method = recode_factor(distance_data_mirna$Method, "hisat2 ns" = "HISAT2", "star enc" = "STAR", "bowtie2 def" = "Bowtie2", "map map" = "vg map", "mpmap_final1 " = "vg mpmap")

distance_data_mirna_old <- distance_data_mirna %>%
  filter(Reads == "ERR187607_uni")

plotDistanceBenchmark(distance_data_mirna_old, wes_cols, "vg_sim_mirna_old.pdf")
  
distance_data_mirna_new <- distance_data_mirna %>%
  filter(Reads == "ENCSR958UOC_rep1_uni")

plotDistanceBenchmark(distance_data_mirna_new, wes_cols, "vg_sim_mirna_new.pdf")

