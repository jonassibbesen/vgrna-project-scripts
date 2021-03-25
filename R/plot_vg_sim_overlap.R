
# plot_vg_sim_overlap.R

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
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  if (grepl("_ovl0_", basename(filename))) {
    
    data <- data %>%
      add_column(Filter = "Unfiltered")
  
  } else if (grepl("_ovl3_", basename(filename))) {
    
    data <- data %>%
      add_column(Filter = "Low quality bases filtered")
    
  } else {
    
    stopifnot(FALSE)
  }
  
  return(data)
}

overlap_threshold <- 90

overlap_data_raw <- map_dfr(list.files(path = "./methods", pattern=".*_ovl.*_vg_.*.txt", full.names = T, recursive = T), parse_file)

overlap_data <- overlap_data_raw %>%
  mutate(Correct = Overlap >= (overlap_threshold / 100)) %>%
  filter(TruthAlignmentLength > 50) 


########


overlap_data_polya <- overlap_data %>%
  filter(Type == "polya_rna")
  
overlap_data_polya$Method <- recode_factor(overlap_data_polya$Method, 
                                     "hisat2" = "HISAT2",
                                     "star" = "STAR",
                                     "map" = "vg map (def)", 
                                     "map_fast" = "vg map", 
                                     "mpmap" = "vg mpmap")

overlap_data_polya <- overlap_data_polya %>%
  filter(Method != "vg map (def)")

overlap_data_polya$FacetCol <- "Simulated reads"
overlap_data_polya$FacetRow <- ""

overlap_data_polya <- overlap_data_polya %>%
  filter(Filter == "Low quality bases filtered")


overlap_data_polya_main <- overlap_data_polya %>%
  filter(Graph != "gencode80") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode80")

overlap_data_polya_main$Graph = recode_factor(overlap_data_polya_main$Graph, 
                                         "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                         "gencode100" = "Spliced reference")

for (reads in unique(overlap_data_polya_main$Reads)) {
  
  overlap_data_polya_main_reads <- overlap_data_polya_main %>%
    filter(Reads == reads)

  plotOverlapBenchmarkMapQ(overlap_data_polya_main_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_main_ovl", overlap_threshold, "_", reads, sep = ""))
}


overlap_data_polya_sj <- overlap_data_polya %>%
  filter(Graph != "gencode100" | Method == "STAR")

overlap_data_polya_sj$Graph = recode_factor(overlap_data_polya_sj$Graph, 
                                       "gencode100" = "All transcripts", 
                                       "1kg_nonCEU_af001_gencode100" = "All transcripts", 
                                       "gencode80" = "80% transcripts", 
                                       "1kg_nonCEU_af001_gencode80" = "80% transcripts")

for (reads in unique(overlap_data_polya_sj$Reads)) {
  
  print(reads)
  
  overlap_data_polya_sj_reads <- overlap_data_polya_sj %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQ(overlap_data_polya_sj_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_sj_ovl", overlap_threshold, "_", reads, sep = ""))
}


# overlap_data_polya_paths <- overlap_data_polya %>%
#   filter(Graph == "1kg_nonCEU_af001_gencode100" | Method == "STAR") 
# 
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "HISAT2",]$Graph <- "Without transcript paths"
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "STAR",]$Graph <- "Without transcript paths"
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "vg map",]$Graph <- "With transcript paths"
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "vg mpmap",]$Graph <- "With transcript paths"
# 
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "map_fast_nopaths",]$Graph <- "Without transcript paths"
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "map_fast_nopaths",]$Method <- "vg map"
# 
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "mpmap_nopaths",]$Graph <- "Without transcript paths"
# overlap_data_polya_paths[overlap_data_polya_paths$Method == "mpmap_nopaths",]$Method <- "vg mpmap"
# 
# overlap_data_polya_paths$Graph <- factor(overlap_data_polya_paths$Graph, levels = c("With transcript paths", "Without transcript paths"))
# plotOverlapBenchmark(overlap_data_polya_paths, wes_cols, "plots/polya_rna/vg_sim_overlap_polya_paths")


########
