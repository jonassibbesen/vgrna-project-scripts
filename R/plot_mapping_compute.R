
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
    add_column(Method = base_split[4]) 
  
  if (base_split[4] == "hisat2" | base_split[4] == "star") {
    
    data <- data %>%
      add_column(Reads = base_split[7]) %>%
      add_column(Graph = base_split[8])
  
  } else {
    
    data <- data %>%
      add_column(Reads = base_split[6]) %>%
      add_column(Graph = base_split[7])   
  }
  
  return(data)    
}

compute_data_hisat2 <- map_dfr(list.files(pattern=".*-hisat2-time-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_star <- map_dfr(list.files(pattern=".*star-time-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_map <- map_dfr(list.files(pattern=".*-map-real-.*log.txt", full.names = T, recursive = T), parse_file)
compute_data_mpmap <- map_dfr(list.files(pattern=".*-mpmap-real-.*gc100-0.*log.txt", full.names = T, recursive = T), parse_file)

compute_data <- bind_rows(compute_data_hisat2, compute_data_star, compute_data_map, compute_data_mpmap) %>%
  add_row(Time = 0, Method = "star", Graph = "na") %>%
  add_row(Time = 0, Method = "star", Graph = "nceu")

compute_data <- compute_data %>%
  filter(Reads == "470")

compute_data$Method = recode_factor(compute_data$Method, "hisat2" = "HISAT2 (incl. SAM->BAM)", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")
compute_data$Graph = recode_factor(compute_data$Graph, "gc100" = "Spliced\nreference", "na" = "Personal\n(NA12878)", "1kg_NA12878_gencode100" = "Personal\n(NA12878)", "nceu" = "1000g\n(no-CEU)")

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

png("real_mapping_compute_srr.png", height = 4, width = 8, units = "in", pointsize = 12, res = 300)
compute_data %>%
  mutate(Time =  Time / 60) %>%
  mutate(Time =  Time / 115359773 * 10 * 10^6) %>%
  ggplot(aes(x = Graph, y = Time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("Minutes per 10M PE reads (16 threads)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))
dev.off()

png("real_mapping_compute_srr_zoom.png", height = 4, width = 8, units = "in", pointsize = 12, res = 300)
compute_data %>%
  mutate(Time =  Time / 60) %>%
  mutate(Time =  Time / 115359773 * 10 * 10^6) %>%
  ggplot(aes(x = Graph, y = Time, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 25), oob = rescale_none) +
  xlab("") +
  ylab("Minutes per 10M PE reads (16 threads)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))
dev.off()
  