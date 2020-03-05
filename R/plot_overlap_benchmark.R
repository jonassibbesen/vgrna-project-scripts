
# plot_overlap_benchmark.R
# Plot mapping overlap benchmark 

rm(list=ls())

library("tidyverse")

args <- commandArgs()
script_dir <- dirname(sub("--file=", "", args[4]))
print(script_dir)

print(args)
system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))

data_dir <- read.csv(args[6], sep = " ", header = F)

setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/null/stats/")

parse_file <- function(file_name) {
  
  fn_split <- strsplit(file_name[1], "_")[[1]]
  
  data <- read_tsv(file_name)
  data <- data %>%
      add_column(Tool = paste(fn_split[2], collapse = "_")) %>%
      add_column(Set = substr(fn_split[length(fn_split)], 1, str_length(fn_split[length(fn_split)]) - 4))
      
  return(data)
}

overlap_data <- map_dfr(list.files(pattern=".*.txt"), parse_file)
overlap_data$Tool = recode_factor(overlap_data$Tool, "hisat2" = "Hisat2", "vg" = "vg mpmap")

overlap_threshold <- 0.9

overlap_data <- overlap_data %>%
  rename(Count = 1) %>%
  mutate(Correct = Overlap >= overlap_threshold) %>% 
  mutate(TP = Count * Correct) %>% 
  mutate(FP = Count * !Correct) 

overlap_data_roc <- overlap_data %>% 
  group_by(Tool, Set, MapQ) %>%
  summarise(TP = sum(TP), FP = sum(FP)) %>% 
  arrange(desc(MapQ), .by_group = T) %>%
  mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
  mutate(N = max(TPcs) + max(FPcs)) %>% 
  mutate(Sensitivity = (FPcs + TPcs) / N, Precision = TPcs / (FPcs + TPcs)) 

pdf(paste(c("roc_overlap", as.integer(overlap_threshold*100), ".pdf"), collapse = ""))
overlap_data_roc %>%
  ggplot(aes(y = Precision, x = Sensitivity, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  facet_wrap(~ Set, ncol=2) +
  xlim(c(0.5, 1)) +
  ylim(c(0.5, 1)) +
  theme_bw() +
  theme(text = element_text(size=16))
dev.off()

overlap_data_mean <- overlap_data %>% 
  mutate(OverlapCount = Overlap * Count) %>%
  group_by(Tool, Set, MapQ) %>%
  summarise(SumCount = sum(Count), SumOverlapCount = sum(OverlapCount)) %>% 
  mutate(MeanOverlap = SumOverlapCount/ SumCount)
  
pdf(paste(c("mean_overlap.pdf"), collapse = ""))
overlap_data_mean %>%
  ggplot(aes(y = MeanOverlap, x = MapQ, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  facet_wrap(~ Set, ncol=2) +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size=16))
dev.off() 
  
overlap_data_mapq <- overlap_data %>% 
  group_by(Tool, Set, MapQ) %>%
  summarise(TP = sum(TP), FP = sum(FP)) %>% 
  mutate(LogMapQError = -1 * MapQ / 10) %>%
  mutate(LogEstError = log10(1 - TP / (FP + TP)))
  
pdf(paste(c("mapq_error.pdf"), collapse = ""))
overlap_data_mapq %>%
  ggplot(aes(y = LogEstError, x = LogMapQError, color = Tool)) +
  geom_line(size = 1) + 
  geom_point(size = 2.5) +
  facet_wrap(~ Set, ncol=2) +
  xlim(c(-6, 0)) +
  ylim(c(-6, 0)) +
  theme_bw() +
  theme(text = element_text(size=16))
dev.off() 
