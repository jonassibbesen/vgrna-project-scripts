
# plot_mapping_compute.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")

#source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r1/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


compute_data <- tibble(Time = numeric(), Method = character(), Reads = character(), Graph = character())

compute_data <- compute_data %>%
  add_row(Time = 45101, Method = "vg map", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 151943, Method = "vg map", Reads = "ENCSR000AED_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 9900, Method = "vg mpmap", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 14598, Method = "vg mpmap", Reads = "ENCSR000AED_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 3040, Method = "HISAT2", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 3110, Method = "HISAT2", Reads = "ENCSR000AED_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 1306, Method = "STAR", Reads = "ENCSR000AED_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = NA, Method = "STAR", Reads = "ENCSR000AED_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 34305, Method = "vg map", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 62988, Method = "vg map", Reads = "CHM13_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 19952, Method = "vg mpmap", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 33694, Method = "vg mpmap", Reads = "CHM13_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 18099, Method = "HISAT2", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = 18192, Method = "HISAT2", Reads = "CHM13_rep1", Graph = "Spliced\npangenome\ngraph") %>%
  add_row(Time = 1963, Method = "STAR", Reads = "CHM13_rep1", Graph = "Spliced\nreference") %>%
  add_row(Time = NA, Method = "STAR", Reads = "CHM13_rep1", Graph = "Spliced\npangenome\ngraph")

parse_ovl_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    summarise(num_reads = sum(Count))
  
  return(data)
}

compute_data$Graph <- factor(compute_data$Graph, levels = c("Spliced\nreference", "Spliced\npangenome\ngraph"))

for (reads in unique(compute_data$Reads)) {
  
  compute_data_reads <- compute_data %>%
    filter(Reads == reads)

  overlap_data <- map_dfr(list.files(path = "./methods", pattern = paste(".*real_r1_", reads, "_exon_ovl_gc.txt.gz", sep = ""), full.names = T, recursive = T), parse_ovl_file) 
  
  compute_data_reads <- compute_data_reads %>%
    mutate(Time = Time * 16) %>%
    mutate(Time = max(overlap_data) / Time / 2) 
  
  compute_data_reads$FacetCol <- "Real reads"
  compute_data_reads$FacetRow <- ""
  
  plotMappingComputeBenchmark(compute_data_reads, wes_cols, paste("plots/real_compute/real_r1_mapping_compute_polya_", reads, sep = ""))
}

########
