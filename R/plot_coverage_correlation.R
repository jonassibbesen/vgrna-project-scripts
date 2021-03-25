
# plot_coverage_correlation.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_types = "iiciii")
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
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
  mutate(MapQ = ifelse(MapQ > 60, 60, MapQ))

coverage_data_pb_mq_corr_list <- list()

for (i in c(0, 1, seq(10, 60, 10))) { 
  
  print(i)
  
  coverage_data_mq <- coverage_data %>%
    mutate(Count = ifelse(MapQ < i, 0, Count)) %>%
    mutate(ReadCoverage = Count * ReadCoverage) %>%
    mutate(BaseCoverage = Count * BaseCoverage) %>%
    group_by(AllelePosition, ExonSize, Type, Reads, Method, Graph) %>%
    summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))
  
  coverage_data_pb_mq <- right_join(pb_coverage, coverage_data_mq, by = c("AllelePosition", "ExonSize")) %>%
    mutate(BaseCoverage.x_norm = BaseCoverage.x / ExonSize) %>%
    mutate(BaseCoverage.y_norm = BaseCoverage.y / ExonSize)
  
  coverage_data_pb_mq_corr_pear <- coverage_data_pb_mq %>%
    group_by(Type, Reads, Method, Graph) %>%
    summarise(num_bases = sum(BaseCoverage.y), Corr = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson")) %>%
    add_column(Threshold = i) %>%
    add_column(cor_type = "Pearson")
    
  coverage_data_pb_mq_corr_list[[as.character(i)]] <- coverage_data_pb_mq_corr_pear
}

coverage_data_pb_mq_corr <- do.call(rbind, coverage_data_pb_mq_corr_list)


########


coverage_data_pb_mq_corr_polya <- coverage_data_pb_mq_corr %>%
  filter(Type == "polya_rna")

coverage_data_pb_mq_corr_polya$Method <- recode_factor(coverage_data_pb_mq_corr_polya$Method, 
                                     "hisat2" = "HISAT2", 
                                     "star" = "STAR", 
                                     "map_fast" = "vg map",
                                     "mpmap" = "vg mpmap")

coverage_data_pb_mq_corr_polya$Graph = recode_factor(coverage_data_pb_mq_corr_polya$Graph, 
                                                   "1kg_nonCEU_af001_gencode100" = "Spliced pangeome graph",
                                                   "gencode100" = "Spliced reference")

coverage_data_pb_mq_corr_polya$FacetCol <- "Real reads"
coverage_data_pb_mq_corr_polya$FacetRow <- ""

for (reads in unique(coverage_data_pb_mq_corr_polya$Reads)) {
  
  coverage_data_pb_mq_corr_polya_reads <- coverage_data_pb_mq_corr_polya %>%
    filter(Reads == reads) %>%
    rename(MapQ = Threshold)
  
  plotIsoSeqCorrelationBenchmark(coverage_data_pb_mq_corr_polya_reads, wes_cols, paste("plots/polya_rna/real_cov_corr_polya_", reads, sep = ""))
}

########
