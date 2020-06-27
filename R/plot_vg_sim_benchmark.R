
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
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping/")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_names = F)
  colnames(data) <- c("Count", "Overlap", "MapQ", "Method")

  data <- data %>%
    select(-Method) %>%
    add_column(Reads = dir_split[6]) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8])

  return(data)
}

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

distance_data_d10 <- map_dfr(list.files(pattern=".*_gam_dist10_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d10$Threshold <- "Distance <= 10"

distance_data_d100 <- map_dfr(list.files(pattern=".*_gam_dist100_.*.txt", full.names = T, recursive = T), parse_file)
distance_data_d100$Threshold <- "Distance <= 100"

distance_data <- bind_rows(distance_data_d10, distance_data_d100) %>%
  mutate(Correct = Overlap > 0.5) 

distance_data$Threshold <- as.factor(distance_data$Threshold)
distance_data$Method = recode_factor(distance_data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")


distance_data_all_sj <- distance_data %>%
  filter(Reads == "SRR1153470_uni") %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") %>%
  filter(Threshold == "Distance <= 10")

distance_data_all_sj$Graph = recode_factor(distance_data_all_sj$Graph, "gencode100" = "Spliced reference", "1kg_NA12878_exons_gencode100" = "Personal (NA12878)", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

plotDistanceBenchmark(distance_data_all_sj, wes_cols, "vg_sim_benchmark_distance_all_sj.pdf")


distance_data_sub_sj <- distance_data %>%
  filter(Reads == "SRR1153470_uni") %>%
  filter(Graph != "gencode100" | Method == "STAR") %>%
  filter(Graph != "1kg_NA12878_exons_gencode100") %>%
  filter(Graph != "1kg_NA12878_gencode100") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") %>%
  filter(Threshold == "Distance <= 10")

distance_data_sub_sj$Graph = recode_factor(distance_data_sub_sj$Graph, "gencode100" = "All splice-junctions", "1kg_nonCEU_af001_gencode100" = "All splice-junctions", "gencode85" = "85% splice-junctions", "1kg_nonCEU_af001_gencode85" = "85% splice-junctions")

plotDistanceBenchmark(distance_data_sub_sj, wes_cols, "vg_sim_benchmark_distance_sub_sj.pdf")


distance_data_gene <- distance_data %>%
  filter(Reads == "SRR1153470_uni") %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "gencode100" | Method == "STAR") %>%
  filter(Graph != "1kg_NA12878_exons_gencode100") %>%
  filter(Graph != "1kg_NA12878_gencode100") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Threshold == "Distance <= 10")

distance_data_gene$Graph = recode_factor(distance_data_gene$Graph, "gencode100" = "Whole genome", "1kg_nonCEU_af001_gencode100" = "Whole genome", "1kg_nonCEU_af001_gencode100_genes" = "Exonic only")

plotDistanceBenchmark(distance_data_gene, wes_cols, "vg_sim_benchmark_distance_gene.pdf")

