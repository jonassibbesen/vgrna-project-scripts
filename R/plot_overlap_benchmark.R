
# plot_overlap_benchmark.R
# Plot mapping overlap benchmark 

rm(list=ls())

library("tidyverse")
library("gridExtra")

# args <- commandArgs()
# script_dir <- dirname(sub("--file=", "", args[4]))
# print(script_dir)
# 
# print(args)
# system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
# system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))
# 
# data_dir <- read.csv(args[6], sep = " ", header = F)

setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/plots/debug/mapping2/1kg_NA12878_exons_gencode/sim_1kg_NA12878_gencode_SRR1153470_avr_rsem//")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_tsv(filename)
  data <- data %>%
    add_column(Tool = dir_split[2]) %>%
    rename(Count = 1)
    
  if ("correct" %in% colnames(data)) {
    
      data <- data %>%
        rename(Overlap = correct) 
  }
  
  if ("mq" %in% colnames(data)) {
    
    data <- data %>%
      rename(MapQ = mq) 
  }
  
  return(data)
}

overlap_data <- map_dfr(list.files(pattern="*.txt", full.names = T, recursive = T), parse_file)
overlap_data$Tool = recode_factor(overlap_data$Tool, "hisat2" = "Hisat2", "vg_mpmap" = "vg mpmap", "vg_mpmap_new_mapq" = "vg mpmap (new mapq)")

overlap_threshold <- 0.9

overlap_data <- overlap_data %>%
  mutate(Correct = Overlap >= overlap_threshold) %>% 
  mutate(TP = Count * Correct) %>% 
  mutate(FP = Count * !Correct) 

overlap_data_roc <- overlap_data %>% 
  group_by(Tool, MapQ) %>%
  summarise(TP = sum(TP), FP = sum(FP)) %>% 
  arrange(desc(MapQ), .by_group = T) %>%
  mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
  mutate(N = max(TPcs) + max(FPcs)) %>% 
  mutate(Sensitivity = (FPcs + TPcs) / N, Precision = TPcs / (FPcs + TPcs)) 

pdf(paste(c("roc_overlap", as.integer(overlap_threshold*100), ".pdf"), collapse = ""))

plot1 <- overlap_data_roc %>%
  ggplot(aes(y = Precision, x = Sensitivity, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  xlim(c(0.75, 1)) +
  ylim(c(0.75, 1)) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(color = guide_legend(nrow = 3)) +
  theme(text = element_text(size=16))

plot2 <- overlap_data_roc %>%
  ggplot(aes(y = log10(1 - Precision), x = Sensitivity, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size=16))

grid.arrange(grobs= lapply(list(plot1, plot2), "+", theme(plot.margin=margin(10,10,10,10))), nrow = 1)

dev.off()

overlap_data_mean <- overlap_data %>% 
  mutate(OverlapCount = Overlap * Count) %>%
  group_by(Tool, MapQ) %>%
  summarise(SumCount = sum(Count), SumOverlapCount = sum(OverlapCount)) %>% 
  mutate(MeanOverlap = SumOverlapCount/ SumCount)
  
pdf(paste(c("mean_overlap.pdf"), collapse = ""))
overlap_data_mean %>%
  ggplot(aes(y = MeanOverlap, x = MapQ, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size=16))
dev.off() 
  
overlap_data_mapq <- overlap_data %>% 
  group_by(Tool, MapQ) %>%
  summarise(TP = sum(TP), FP = sum(FP)) %>% 
  mutate(LogMapQError = -1 * MapQ / 10) %>%
  mutate(LogEstError = log10(1 - TP / (FP + TP)))
  
pdf(paste(c("mapq_error", as.integer(overlap_threshold*100), ".pdf"), collapse = ""))
overlap_data_mapq %>%
  ggplot(aes(y = LogEstError, x = LogMapQError, color = Tool)) +
  geom_abline(intercept = 0) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  xlim(c(-6, 0)) +
  ylim(c(-6, 0)) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=16))
dev.off() 

overlap_data_frac <- overlap_data %>%
  mutate(Correct = Overlap == 1) %>% 
  mutate(TP = Count * Correct) %>% 
  mutate(FP = Count * !Correct) %>%
  group_by(Tool, MapQ) %>%
  summarise(TP = sum(TP), FP = sum(FP)) %>% 
  mutate(FracPerfect = TP / (FP + TP))

pdf(paste(c("frac_perfect.pdf"), collapse = ""))
overlap_data_frac %>%
  ggplot(aes(y = FracPerfect, x = MapQ, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size=16))
dev.off() 
