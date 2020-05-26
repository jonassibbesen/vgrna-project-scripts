
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
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/plots/mapping/")

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
  time = 0
  
  if (dir_split[2] == "hisat2" | dir_split[2] == "star") {
    
    if (dir_split[4] == "real_SRR1153470") {
      
      time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1])
      
    } else if (dir_split[4] == "sim_SRR1153470_uni_rsem") {
      
      time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1]) + convertTimeLine(grep('Elapsed', readLines(filename), value = T)[3])
      
    } else if (dir_split[4] == "sim_SRR1153470_uni_vg") {
      
      time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[2]) + convertTimeLine(grep('Elapsed', readLines(filename), value = T)[5])
    }
    
  } else {
    
    if (dir_split[4] == "real_SRR1153470") {
      
      time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1])
      
    } else if (dir_split[4] == "sim_SRR1153470_uni_rsem") {  
      
      time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1]) + convertTimeLine(grep('Elapsed', readLines(filename), value = T)[2])
      
    } else if (dir_split[4] == "sim_SRR1153470_uni_vg") {
      
      time <- convertTimeLine(grep('Elapsed', readLines(filename), value = T)[1]) + convertTimeLine(grep('Elapsed', readLines(filename), value = T)[4])
    }     
  }
  
  data <- data_frame(Time = time)
  
  if (dir_split[2] == "vg") {
    
    data <- data %>%
      add_column(Reads = dir_split[4]) %>%
      add_column(Method = dir_split[5])
    
  } else {
    
    data <- data %>%
      add_column(Method = dir_split[2]) %>%
      add_column(Reads = dir_split[4])
  }
  
  if (grepl("nonceu", basename(filename), fixed = TRUE)) {
    
    data <- data %>%
      add_column(Graph = "1000g")    
  
  } else if (grepl("nceu", basename(filename), fixed = TRUE)) {
      
      data <- data %>%
        add_column(Graph = "1000g")    
      
  } else if (grepl("na12878", basename(filename), fixed = TRUE)) {
    
    data <- data %>%
      add_column(Graph = "NA12878")    
  
  } else {
    
    data <- data %>%
      add_column(Graph = "Linear")       
  }

  data$Method = recode_factor(data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap", "mpmap_old" = "vg mpmap\n(old)", "mpmap_fast" = "vg mpmap\n(fast)", "mpmap_fast2" = "vg mpmap\n(fast2)", "mpmap_fast3" = "vg mpmap\n(fast3)", "mpmap_fast4" = "vg mpmap\n(fast4)", "mpmap_fast5" = "vg mpmap\n(fast5)", "mpmap_fast6" = "vg mpmap\n(fast6)", "mpmap_fast7" = "vg mpmap\n(fast7)")
  data$Reads = recode_factor(data$Reads, "real_SRR1153470" = "Real", "sim_SRR1153470_uni_rsem" = "Sim (RSEM)", "sim_SRR1153470_uni_vg" = "Sim (vg)")
  
  return(data)    
}

compute_data <- map_dfr(list.files(pattern=".*reads.*gc100.*log.txt", full.names = T, recursive = T), parse_file)

compute_data <- compute_data %>%
  filter(Graph == "1000g") %>%
  filter(Reads != "sim_SRR1153470_rsem")  %>%
  filter(Reads == "Real")  %>%
  mutate(Time = Time / 60)

compute_data <- compute_data %>%
  spread(Method, Time) %>%
  replace_na(list(`vg mpmap\n(old)` = 0)) %>%
  gather(HISAT2, `vg map`, `vg mpmap`, `vg mpmap\n(old)`, `vg mpmap\n(fast)`, `vg mpmap\n(fast2)`, `vg mpmap\n(fast3)`, `vg mpmap\n(fast4)`, `vg mpmap\n(fast6)`, `vg mpmap\n(fast7)`, key = "Method", value = "Time")

wes_cols <- wes_palette("Darjeeling1")[c(3,5)]

pdf("compute_times.pdf", width = 10)

compute_data %>% 
  ggplot(aes(y = Time, x = Method, fill = Graph)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(rows = vars(Reads), scales = "free") + 
  scale_fill_manual(values = wes_cols) +
  ylab("Compute time (minutes)") +
  xlab("") +
  theme_bw() +
  theme(text = element_text(size=16))

compute_data <- compute_data %>%
  spread(Method, Time) %>%
  mutate(`vg map` = `vg map` / HISAT2) %>%
  mutate(`vg mpmap` = `vg mpmap` / HISAT2) %>%
  mutate(`vg mpmap\n(old)` = `vg mpmap\n(old)` / HISAT2) %>%
  mutate(`vg mpmap\n(fast)` = `vg mpmap\n(fast)` / HISAT2) %>%
  mutate(`vg mpmap\n(fast2)` = `vg mpmap\n(fast2)` / HISAT2) %>%
  mutate(`vg mpmap\n(fast3)` = `vg mpmap\n(fast3)` / HISAT2) %>%
  mutate(`vg mpmap\n(fast4)` = `vg mpmap\n(fast4)` / HISAT2) %>%
  mutate(`vg mpmap\n(fast6)` = `vg mpmap\n(fast6)` / HISAT2) %>%
  mutate(`vg mpmap\n(fast7)` = `vg mpmap\n(fast7)` / HISAT2) %>%
  mutate(HISAT2 = HISAT2 / HISAT2) %>%
  gather(HISAT2, `vg map`, `vg mpmap`, `vg mpmap\n(old)`, `vg mpmap\n(fast)`, `vg mpmap\n(fast2)`, `vg mpmap\n(fast3)`, `vg mpmap\n(fast4)`, `vg mpmap\n(fast6)`, `vg mpmap\n(fast7)`, key = "Method", value = "Time")

compute_data %>% 
  ggplot(aes(y = Time, x = Method, fill = Graph)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(rows = vars(Reads)) + 
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 15),oob = rescale_none) +
  ylab("Times slower than HISAT2") +
  xlab("") +
  theme_bw() +
  theme(text = element_text(size=16))

dev.off()
