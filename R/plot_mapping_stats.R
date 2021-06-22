
# plot_mapping_stats.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])

  return(data)
}

mapping_data <- map_dfr(list.files(path = "./methods", pattern=".*_exon_ovl_gc.*.txt", full.names = T, recursive = T), parse_file)


########


mapping_data_polya <- mapping_data %>%
  filter(Type == "polya_rna")

mapping_data_polya$Method <- recode_factor(mapping_data_polya$Method, 
                                           "hisat2" = "HISAT2",
                                           "star" = "STAR",
                                           "map_fast" = "vg map", 
                                           "mpmap" = "vg mpmap")

mapping_data_polya$Graph = recode_factor(mapping_data_polya$Graph, 
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                         "1kg_all_af001_gencode100" = "Spliced pangenome graph")

mapping_data_polya_stats <- mapping_data_polya %>%
  mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
  mutate(MapQ0 = Count * (MapQ >= 0)) %>% 
  mutate(MapQ1 = Count * (MapQ >= 1)) %>% 
  group_by(Reads, Method, Graph) %>%
  summarise(count = sum(Count), MapQ0 = sum(MapQ0), MapQ1 = sum(MapQ1)) %>%
  mutate(MapQ0_frac = MapQ0 / count) %>%
  mutate(MapQ1_frac = MapQ1 / count) %>%
  gather("MapQ0_frac", "MapQ1_frac", key = "Filter", value = "Frac")

for (reads in unique(mapping_data_polya_stats$Reads)) {
  
  mapping_data_polya_stats_reads <- mapping_data_polya_stats %>%
    filter(Reads == reads)
  
  mapping_data_polya_stats_reads <- mapping_data_polya_stats_reads %>%
    ungroup() %>%
    add_row(Reads = reads, Method = "STAR", Graph = "Spliced pangenome graph", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
    add_row(Reads = reads, Method = "STAR", Graph = "Spliced pangenome graph", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0)

    
  mapping_data_polya_stats_reads$Graph = recode_factor(mapping_data_polya_stats_reads$Graph, 
                                           "Spliced reference" = "Spliced\nreference",
                                           "Spliced pangenome graph" = "Spliced pan-\ngenome graph")
  
  mapping_data_polya_stats_reads$Filter <- recode_factor(mapping_data_polya_stats_reads$Filter, 
                                                   "MapQ0_frac" = "All",
                                                   "MapQ1_frac" = "MapQ > 0")
  
  mapping_data_polya_stats_reads$FacetCol <- "Real reads"
  mapping_data_polya_stats_reads$FacetRow <- ""
  
  plotMappingStatsBenchmark(mapping_data_polya_stats_reads, wes_cols, paste("plots/polya_rna/real_stats_polya_bar_", reads, sep = ""))
}

########
