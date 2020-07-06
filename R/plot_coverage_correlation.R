
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
    add_column(Reads = dir_split[6]) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8])
  
  return(data)
}

pb_coverage <- read_table2("alignments/ENCSR706ANY/ENCSR706ANY_exon_cov_bam.txt", col_types = "iiciii")
pb_coverage <- pb_coverage %>%
  mutate(Count = ifelse(MapQ < 30, 0, Count)) %>%
  mutate(ReadCoverage = Count * ReadCoverage) %>%
  mutate(BaseCoverage = Count * BaseCoverage) %>%
  group_by(AllelePosition, ExonSize) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

coverage_data <- map_dfr(list.files(pattern=".*_real_.*_exon_cov_ENCSR706ANY.txt", full.names = T, recursive = T), parse_file)

n_reads = 230719546

coverage_data_pb_mq_cor_list <- list()

for (i in seq(0, 60, 1)) { 
  
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
    summarise(sens = sum(BaseCoverage.y) / (n_reads * 101), cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Pearson")
    
  coverage_data_pb_mq_cor_spea <- coverage_data_pb_mq %>%
    group_by(Reads, Method, Graph) %>%
    summarise(sens = sum(BaseCoverage.y) / (n_reads * 101), cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "spearman")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Spearman")
    
    coverage_data_pb_mq_cor_list[[as.character(i)]] <- rbind(coverage_data_pb_mq_cor_pear, coverage_data_pb_mq_cor_spea)
}

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

coverage_data_pb_mq_cor_data <- do.call(rbind, coverage_data_pb_mq_cor_list)

coverage_data_pb_mq_cor_data$Method = recode_factor(coverage_data_pb_mq_cor_data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")
coverage_data_pb_mq_cor_data$Graph = recode_factor(coverage_data_pb_mq_cor_data$Graph, "gencode100" = "Spliced reference", "1kg_NA12878_exons_gencode100" = "Personal (NA12878)", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

pdf("real_benchmark_exon_cor_srr.pdf", width = 10)
coverage_data_pb_mq_cor_data %>%
  ggplot(aes(y = cor, x = sens, color = Method, linetype = Graph, shape = Graph)) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) +
  scale_color_manual(values = wes_cols) +
  facet_grid(cols = vars(cor_type)) +
  xlab("Fraction mapped bases overlapping Iso-Seq reads") +
  ylab("Iso-Seq coverage correlation (exon average)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()


