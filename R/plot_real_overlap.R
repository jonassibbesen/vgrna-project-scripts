
# plot_real_overlap.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping/")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  return(data)
}

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

overlap_data_o1 <- map_dfr(list.files(pattern=".*_real_.*_ol_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o1$Threshold <- "Overlap >= 1%"

overlap_data_o1 <- overlap_data_o1 %>%
  mutate(Correct = Overlap > 0.01)

overlap_data_o90 <- map_dfr(list.files(pattern=".*_real_.*_ol_.*.txt", full.names = T, recursive = T), parse_file)
overlap_data_o90$Threshold <- "Overlap >= 90%"

overlap_data_o90 <- overlap_data_o90 %>%
  mutate(Correct = Overlap > 0.90)

overlap_data <- bind_rows(overlap_data_o90, overlap_data_o1) 

overlap_data$Threshold <- factor(overlap_data$Threshold, levels = c("Overlap >= 90%", "Overlap >= 1%"))
overlap_data$Method = recode_factor(overlap_data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")

overlap_data$Graph = recode_factor(overlap_data$Graph, "gencode100" = "Spliced reference", "1kg_NA12878_exons_gencode100" = "Personal (NA12878)", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

overlap_data <- overlap_data %>%
  filter(Reads != "ENCSR000AED_rep2")

overlap_data$Reads = recode_factor(overlap_data$Reads, "SRR1153470" = "Training set", "ENCSR000AED_rep1" = "Test set")


plotOverlapBenchmark(overlap_data, wes_cols, "real_benchmark_overlap")

