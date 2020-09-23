
# plot_mapping_memory.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping/")

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
    add_column(Reads = base_split[6]) %>%
    add_column(Graph = base_split[7])
  
  if (grepl("-gs-", basename(filename))) {
    
    data <- data %>%
      mutate(Graph = paste(Graph, "gs", sep = "_"))
  
  } else if (grepl("-gt10-", basename(filename))) {
  
    data <- data %>%
      mutate(Graph = paste(Graph, "gt10", sep = "_"))
  }
  
  if (grepl("-f-", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "fast", sep = "_")) %>%
      mutate(Reads = base_split[7]) %>%
      mutate(Graph = base_split[8])
  }
  
  return(data)    
}

memory_data_hisat2 <- map_dfr(list.files(pattern="-hisat2-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_star <- map_dfr(list.files(pattern=".*star-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_map <- map_dfr(list.files(pattern=".*-map-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_map_fast <- map_dfr(list.files(pattern=".*-map-f-real-.*log.txt", full.names = T, recursive = T), parse_file)
memory_data_mpmap <- map_dfr(list.files(pattern=".*-mpmap-real-.*log.txt", full.names = T, recursive = T), parse_file)

memory_data <- bind_rows(memory_data_hisat2, memory_data_star, memory_data_map, memory_data_map_fast, memory_data_mpmap) %>%
  filter(Reads == data_set2) %>%
  filter(Graph != "nceu_gs") %>%
  filter(Graph != "nceu_gt10")

memory_data <- memory_data %>%
  add_row(Mem = 0, Method = "star", Graph = "nceu") 

memory_data$Method = recode_factor(memory_data$Method, 
                                    "hisat2" = "HISAT2", 
                                    "star" = "STAR", 
                                    "map" = "vg map (def)", 
                                    "map_fast" = "vg map", 
                                    "mpmap" = "vg mpmap")

memory_data <- memory_data %>%
  filter(Method != "vg map (def)")

memory_data$Graph = recode_factor(memory_data$Graph, 
                                   "gc100" = "Spliced\nreference",
                                   "nceu" = "1000g\n(no-CEU)")

pdf("plots/polya_rna/real_mapping_memory.pdf", height = 4, width = 4, pointsize = 12)
memory_data %>%
  ggplot(aes(x = Graph, y = Mem, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 60), oob = rescale_none) +
  xlab("") +
  ylab("Maximum memory usage (GiB)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=12))
dev.off()
