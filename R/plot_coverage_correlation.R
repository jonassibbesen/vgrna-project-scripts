
# plot_coverage_correlation.R

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
  
  data <- read_table2(filename, col_types = "iiciii")
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])

  return(data)
}

pb_coverage <- read_table2("alignments/ENCSR706ANY/ENCSR706ANY_mq30_exon_cov_bam.txt.gz", col_types = "iiciii")
pb_coverage <- pb_coverage %>%
  mutate(ReadCoverage = Count * ReadCoverage) %>%
  mutate(BaseCoverage = Count * BaseCoverage) %>%
  group_by(AllelePosition, ExonSize) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_real_.*cov_ENCSR706ANY_mq30.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Reads == data_set3) 

coverage_data_pb_mq_cor_list <- list()

for (i in c(0, 1, seq(10, 60, 10), 255)) { 
  
  print(i)
  
  coverage_data_mq <- coverage_data %>%
    mutate(Count = ifelse(MapQ < i, 0, Count)) %>%
    mutate(ReadCoverage = Count * ReadCoverage) %>%
    mutate(BaseCoverage = Count * BaseCoverage) %>%
    group_by(AllelePosition, ExonSize, Reads, Method, Graph) %>%
    summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))
  
  coverage_data_pb_mq <- right_join(pb_coverage, coverage_data_mq, by = c("AllelePosition", "ExonSize")) %>%
    mutate(BaseCoverage.x_norm = BaseCoverage.x / ExonSize) %>%
    mutate(BaseCoverage.y_norm = BaseCoverage.y / ExonSize)
  
  coverage_data_pb_mq_cor_pear <- coverage_data_pb_mq %>%
    group_by(Reads, Method, Graph) %>%
    summarise(num_bases = sum(BaseCoverage.y), cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Pearson")
    
    coverage_data_pb_mq_cor_list[[as.character(i)]] <- coverage_data_pb_mq_cor_pear
}

parse_ovl_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    mutate(map_length = Length - InsertionLength - SoftClipLength) %>%
    group_by(Reads, Method, Graph, MapQ) %>%
    summarise(num_reads = sum(Count), num_read_bases = sum(Length * Count), num_mapped_bases = sum(IsMapped * Count * map_length)) %>%
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(num_reads_cs = cumsum(num_reads), num_read_bases_cs = cumsum(num_read_bases), num_mapped_bases_cs = cumsum(num_mapped_bases))
    
  return(data)
}

overlap_data <- map_dfr(list.files(path = "./methods", pattern=".*_exon_ovl_ENCSR706ANY_mq30.*.txt", full.names = T, recursive = T), parse_ovl_file) %>%
  rename(Threshold = MapQ)

coverage_data_pb_mq_cor_data <- do.call(rbind, coverage_data_pb_mq_cor_list) %>%
  filter(!is.na(cor)) %>%
  left_join(overlap_data, by = c("Reads", "Method", "Graph", "Threshold")) %>%
  group_by(Reads, Method, Graph, cor_type) %>%
  arrange(desc(Threshold), .by_group = T) %>%
  filter(!is.na(num_reads)) %>%
  mutate(frac_bases = num_bases / max(num_read_bases_cs))

coverage_data_pb_mq_cor_data$Method <- recode_factor(coverage_data_pb_mq_cor_data$Method, 
                                     "hisat2" = "HISAT2", 
                                     "star" = "STAR", 
                                     "map_fast" = "vg map",
                                     "mpmap" = "vg mpmap")

coverage_data_pb_mq_cor_data$Graph = recode_factor(coverage_data_pb_mq_cor_data$Graph, 
                                   "gencode100" = "Spliced reference",
                                   "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

coverage_data_pb_mq_cor_data <- coverage_data_pb_mq_cor_data %>%
  mutate(Threshold = ifelse(Threshold > 60, 60, Threshold))

pdf("plots/polya_rna/real_cov_correlation.pdf", height = 4, width = 6, pointsize = 12)
coverage_data_pb_mq_cor_data %>%
  ggplot(aes(y = cor, x = Threshold, color = Method, linetype = Graph, shape = Graph, label = Threshold)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  #geom_point(data = subset(coverage_data_pb_mq_cor_data, Threshold == 0 | Threshold == 1 | Threshold == 60 | Threshold == 255), size = 2.5, alpha = 1) +
  #geom_text_repel(data = subset(coverage_data_pb_mq_cor_data, Threshold == 0 | Threshold == 1 | Threshold == 60 | Threshold == 255), size = 3.5, fontface = 2) + 
  scale_color_manual(values = wes_cols) +
  xlab("Mapping quality threshold") +
  ylab("Iso-Seq exon coverage correlation") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))
dev.off()

