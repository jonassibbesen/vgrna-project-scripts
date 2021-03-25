
# plot_mapping_memory.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)


convertMemoryLine <- function(memory_line) {

  return(as.double(strsplit(memory_line, " ")[[1]][6]) / 10^6)
}

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  base_split <- strsplit(basename(filename), "-")[[1]]
  
  mem <- convertMemoryLine(grep('Maximum', readLines(filename), value = T)[1])
  
  data <- data_frame(Mem = mem)
  data <- data %>%
    add_column(Method = base_split[4]) %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = base_split[6]) %>%
    add_column(Graph = base_split[7])
  
  if (grepl("-f-", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "fast", sep = "_")) %>%
      mutate(Reads = base_split[7]) %>%
      mutate(Graph = base_split[8])
  }
  
  if (grepl("-gs-", basename(filename))) {
    
    data <- data %>%
      mutate(Graph = paste(Graph, "gs", sep = "_"))
  
  } else if (grepl("-gt10-", basename(filename))) {
  
    data <- data %>%
      mutate(Graph = paste(Graph, "gt10", sep = "_"))
  }
  
  return(data)    
}

memory_data_hisat2 <- map_dfr(list.files(pattern="-hisat2-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_star <- map_dfr(list.files(pattern=".*star-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_map_fast <- map_dfr(list.files(pattern=".*-map-f-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_mpmap <- map_dfr(list.files(pattern=".*-mpmap-real-.*log.txt", full.names = T, recursive = T), parse_file)

memory_data <- bind_rows(memory_data_hisat2, memory_data_star, memory_data_map_fast, memory_data_mpmap)


########


memory_data_polya <- memory_data %>%
  filter(Type == "polya_rna")

memory_data_polya$Method = recode_factor(memory_data_polya$Method, 
                                    "hisat2" = "HISAT2", 
                                    "star" = "STAR", 
                                    "map_fast" = "vg map", 
                                    "mpmap" = "vg mpmap")

memory_data_polya$Reads = recode_factor(memory_data_polya$Reads, 
                                         "aed1" = "ENCSR000AED_rep1", 
                                         "t2t1" = "CHM13_rep1", 
                                         "470" = "SRR1153470")

for (reads in unique(memory_data_polya$Reads)) {
  
  memory_data_polya_reads <- memory_data_polya %>%
    filter(Reads == reads)
  
  if (reads == "ENCSR000AED_rep1" | reads == "SRR1153470") {
    
    memory_data_polya_reads <- memory_data_polya_reads %>%
      filter(Graph != "all") %>%
      add_row(Mem = 0, Method = "STAR", Graph = "nceu")
  }
  
  if (reads == "CHM13_rep1") {
    
    memory_data_polya_reads <- memory_data_polya_reads %>%
      add_row(Mem = 0, Method = "STAR", Graph = "all")
  }
  
  memory_data_polya_reads$Graph = recode_factor(memory_data_polya_reads$Graph, 
                                                 "gc100" = "Spliced\nreference",
                                                 "nceu" = "Spliced pan-\ngenome graph",
                                                 "all" = "Spliced pan-\ngenome graph")
  
  memory_data_polya_reads$FacetCol <- "Real reads"
  memory_data_polya_reads$FacetRow <- ""
  
  plotMappingMemoryBenchmark(memory_data_polya_reads, wes_cols, paste("plots/polya_rna/real_mapping_memory_polya_", reads, sep = ""))
}

########
