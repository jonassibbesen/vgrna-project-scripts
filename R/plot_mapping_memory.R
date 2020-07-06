
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
    add_column(Method = base_split[4]) 
  
  if (base_split[4] == "hisat2" | base_split[4] == "star") {
    
    data <- data %>%
      add_column(Graph = base_split[8])
    
  } else {
    
    data <- data %>%
      add_column(Graph = base_split[7])   
  }
  
  return(data)    
}


compute_data_hisat2 <- map_dfr(list.files(pattern=".*-hisat2-bam-real-470.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_star <- map_dfr(list.files(pattern=".*star-bam-real-470.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_map <- map_dfr(list.files(pattern=".*-map-real-470.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_mpmap <- map_dfr(list.files(pattern=".*-mpmap-real-470.*gc100-0.*log.txt", full.names = T, recursive = T), parse_file)

compute_data <- bind_rows(compute_data_hisat2, compute_data_star, compute_data_map, compute_data_mpmap) %>%
  add_row(Mem = 0, Method = "star", Graph = "na") %>%
  add_row(Mem = 0, Method = "star", Graph = "nceu")

compute_data$Method = recode_factor(compute_data$Method, "hisat2" = "HISAT2 (incl. SAM->BAM)", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")
compute_data$Graph = recode_factor(compute_data$Graph, "gc100" = "Spliced\nreference", "na" = "Personal\n(NA12878)", "1kg_NA12878_gencode100" = "Personal\n(NA12878)", "nceu" = "1000g\n(no-CEU)")

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

pdf("real_mapping_memory_srr.pdf", width = 8)
compute_data %>%
  ggplot(aes(x = Graph, y = Mem, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 100), oob = rescale_none) +
  xlab("") +
  ylab("Maximum memory usage (GiB)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()
