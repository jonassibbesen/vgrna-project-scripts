
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

coverage_data <- map_dfr(list.files(pattern=".*_real_.*cov_ENCSR706ANY_mq30.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Reads == data_set3) 

coverage_data_pb_mq_cor_list <- list()

for (i in c(0, 1, seq(5, 60, 5), 255)) { 
  
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
    summarise(sens = sum(BaseCoverage.y) / 101, cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Pearson")
   
  coverage_data_pb_mq_cor_log_pear <- coverage_data_pb_mq %>%
    group_by(Reads, Method, Graph) %>%
    summarise(sens = sum(BaseCoverage.y) / 101, cor = cor(log(BaseCoverage.x_norm + 1), log(BaseCoverage.y_norm + 1), method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "LogPearson")
   
  coverage_data_pb_mq_cor_spea <- coverage_data_pb_mq %>%
    group_by(Reads, Method, Graph) %>%
    summarise(sens = sum(BaseCoverage.y) / 101, cor = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "spearman")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Spearman")
  
  coverage_data_pb_mq_cor_exp <- coverage_data_pb_mq %>%
    group_by(Reads, Method, Graph) %>%
    summarise(sens = sum(BaseCoverage.y) / 101, cor = mean(BaseCoverage.y > 0)) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Expressed")
    
    coverage_data_pb_mq_cor_list[[as.character(i)]] <- rbind(coverage_data_pb_mq_cor_pear, coverage_data_pb_mq_cor_log_pear, coverage_data_pb_mq_cor_spea, coverage_data_pb_mq_cor_exp)
}

coverage_data_pb_mq_cor_data <- do.call(rbind, coverage_data_pb_mq_cor_list) %>%
  filter(!is.na(cor))

coverage_data_pb_mq_cor_data[coverage_data_pb_mq_cor_data$Reads == data_set3,]$sens <- coverage_data_pb_mq_cor_data[coverage_data_pb_mq_cor_data$Reads == data_set3,]$sens / (2 *  num_reads[[data_set3]])

coverage_data_pb_mq_cor_data$Method <- recode_factor(coverage_data_pb_mq_cor_data$Method, 
                                     "hisat2" = "HISAT2", 
                                     "star" = "STAR", 
                                     "map" = "vg map (def)", 
                                     "map_fast" = "vg map",
                                     "mpmap" = "vg mpmap")

coverage_data_pb_mq_cor_data <- coverage_data_pb_mq_cor_data %>%
  filter(Method != "vg map (def)")

coverage_data_pb_mq_cor_data$Graph = recode_factor(coverage_data_pb_mq_cor_data$Graph, 
                                   "gencode100" = "Spliced reference",
                                   "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)",
                                   "1kg_nonCEU_af001_gencode100_gtex10s2r8e1g" = "1000g (GTEx)")

coverage_data_pb_mq_cor_data_type <- coverage_data_pb_mq_cor_data %>%
  filter(cor_type != "LogPearson") %>%
  filter(cor_type != "Expressed") %>%
  filter(cor_type != "Spearman")

pdf("plots/polya_rna/real_cov_correlation.pdf", height = 5, width = 7, pointsize = 12)
coverage_data_pb_mq_cor_data_type %>%
  ggplot(aes(y = cor, x = sens, color = Method, linetype = Graph, shape = Graph, label = Threshold)) +
  geom_line(size = 1) +
  geom_point(data = subset(coverage_data_pb_mq_cor_data_type, Threshold == 0 | Threshold == 1 | Threshold == 60 | Threshold == 255), size = 2.75, alpha = 1) +
  geom_text_repel(data = subset(coverage_data_pb_mq_cor_data_type, Threshold == 0 | Threshold == 1 | Threshold == 60 | Threshold == 255), size = 4, fontface = 2) + 
  scale_color_manual(values = wes_cols) +
  xlab("Fraction mapped bases overlapping Iso-Seq reads") +
  ylab("Iso-Seq exon coverage correlation") +
  guides(color = FALSE) +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=16))
dev.off()

