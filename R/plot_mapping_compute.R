
# plot_mapping_compute.R

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

convertTimeLine <- function(time_line) {
  
  time <- 0
  time_idx = 0
  
  for (value in rev(strsplit(strsplit(time_line, " ")[[1]][8], ":")[[1]])) {
    
    time = time + as.double(value) * 60^time_idx
    time_idx = time_idx + 1
  }
  
  return(time)
}

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  base_split <- strsplit(basename(filename), "-")[[1]]
  
  time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1])

  if (base_split[4] == "hisat2") {
    
    time <- time + convertTimeLine(grep('Elapsed', readLines(filename), value = T)[2])
  }
  
  data <- data_frame(Time = time)
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

compute_data_hisat2 <- map_dfr(list.files(pattern="-hisat2-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_star <- map_dfr(list.files(pattern=".*star-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_map <- map_dfr(list.files(pattern=".*-map-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_map_fast <- map_dfr(list.files(pattern=".*-map-f-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_mpmap <- map_dfr(list.files(pattern=".*-mpmap-real-.*log.txt", full.names = T, recursive = T), parse_file)

compute_data <- bind_rows(compute_data_hisat2, compute_data_star, compute_data_map, compute_data_map_fast, compute_data_mpmap) %>%
   filter(Reads == data_set2) %>%
   filter(Graph != "nceu_gs") %>%
   filter(Graph != "nceu_gt10")
  
compute_data <- compute_data %>%
  add_row(Time = 0, Method = "star", Graph = "nceu") 

compute_data$Method = recode_factor(compute_data$Method, 
                                   "hisat2" = "HISAT2", 
                                   "star" = "STAR", 
                                   "map" = "vg map (def)", 
                                   "map_fast" = "vg map", 
                                   "mpmap" = "vg mpmap")

compute_data <- compute_data %>%
  filter(Method != "vg map (def)")

compute_data$Graph = recode_factor(compute_data$Graph, 
                                  "gc100" = "Spliced\nreference",
                                  "nceu" = "1000g\n(no-CEU)")

pdf("plots/polya_rna/real_mapping_compute.pdf", height = 4, width = 4, pointsize = 12)
compute_data %>%
  mutate(Time =  Time / 60) %>%
  mutate(Time =  Time / num_reads[[data_set2]] * 10 * 10^6) %>%
  ggplot(aes(x = Graph, y = Time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("Mapping time per 10M reads (16 threads)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=12))
dev.off()

pdf("plots/polya_rna/real_mapping_compute_zoom.pdf", height = 4, width = 4, pointsize = 12)
compute_data %>%
  mutate(Time =  Time / 60) %>%
  mutate(Time =  Time / num_reads[[data_set2]] * 10 * 10^6) %>%
  ggplot(aes(x = Graph, y = Time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 50), oob = rescale_none) +
  xlab("") +
  ylab("Mapping time per 10M reads (16 threads)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=12))
dev.off()
  