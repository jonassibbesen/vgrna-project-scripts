
# plot_vg_sim_overlap.R

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
  
  data <- read_table2(filename)
  data <- data %>%
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

overlap_data_raw <- map_dfr(list.files(pattern=".*_ovl.*_vg_.*.txt", full.names = T, recursive = T), parse_file)

overlap_data_o10 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.1) %>%
  add_column(Threshold = "Overlap >= 10%") 

overlap_data_o90 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.9) %>%
  add_column(Threshold = "Overlap >= 90%") 

overlap_data <- bind_rows(overlap_data_o10, overlap_data_o90) %>%
  filter(TruthAlignmentLength >= 10) 


########


overlap_data_polya <- overlap_data %>%
  filter(Reads == data_set1) 

overlap_data_polya$Method <- recode_factor(overlap_data_polya$Method, 
                                     "hisat2" = "HISAT2",
                                     "star" = "STAR",
                                     "map" = "vg map (def)", 
                                     "map_fast" = "vg map", 
                                     "mpmap" = "vg mpmap")

overlap_data_polya <- overlap_data_polya %>%
  filter(Method != "vg map (def)")

overlap_data_polya$Reads <- recode_factor(overlap_data_polya$Reads, 
                                     "SRR1153470_uni" = "Training set", 
                                     "ENCSR000AED_rep1_uni" = "Test set")
  
overlap_data_polya$FacetCol <- recode_factor(overlap_data_polya$Threshold, 
                                             "Overlap >= 90%" = "Overlap >= 90%", 
                                             "Overlap >= 10%" = "Overlap >= 10%")

overlap_data_polya$FacetRow <- recode_factor(overlap_data_polya$Filter, 
                                             "Unfiltered" = "Unfiltered", 
                                             "Low quality bases filtered (< 4)" = "Low quality bases filtered (< 4)")


overlap_data_polya_main <- overlap_data_polya %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") %>%
  filter(Method != "map_fast_nopaths") %>%
  filter(Method != "mpmap_nopaths") %>%
  filter(Method != "mpmap_nopaths_gamp") 

overlap_data_polya_main$Graph = recode_factor(overlap_data_polya_main$Graph, 
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

plotOverlapBenchmark(overlap_data_polya_main, wes_cols, "plots/polya_rna/vg_sim_overlap_polya_main")


overlap_data_polya_sj <- overlap_data_polya %>%
  filter(Graph != "gencode100" | Method == "STAR") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") %>%
  filter(Method != "map_fast_nopaths") %>%
  filter(Method != "mpmap_nopaths") %>%
  filter(Method != "mpmap_nopaths_gamp") 

overlap_data_polya_sj$Graph = recode_factor(overlap_data_polya_sj$Graph, 
                                       "gencode100" = "All transcripts", 
                                       "1kg_nonCEU_af001_gencode100" = "All transcripts", 
                                       "gencode85" = "85% transcripts", 
                                       "1kg_nonCEU_af001_gencode85" = "85% transcripts")

plotOverlapBenchmark(overlap_data_polya_sj, wes_cols, "plots/polya_rna/vg_sim_overlap_polya_sj")


overlap_data_polya_paths <- overlap_data_polya %>%
  filter(Graph == "1kg_nonCEU_af001_gencode100" | Method == "STAR") 

overlap_data_polya_paths[overlap_data_polya_paths$Method == "HISAT2",]$Graph <- "Without transcript paths"
overlap_data_polya_paths[overlap_data_polya_paths$Method == "STAR",]$Graph <- "Without transcript paths"
overlap_data_polya_paths[overlap_data_polya_paths$Method == "vg map",]$Graph <- "With transcript paths"
overlap_data_polya_paths[overlap_data_polya_paths$Method == "vg mpmap",]$Graph <- "With transcript paths"

overlap_data_polya_paths[overlap_data_polya_paths$Method == "map_fast_nopaths",]$Graph <- "Without transcript paths"
overlap_data_polya_paths[overlap_data_polya_paths$Method == "map_fast_nopaths",]$Method <- "vg map"

overlap_data_polya_paths[overlap_data_polya_paths$Method == "mpmap_nopaths",]$Graph <- "Without transcript paths"
overlap_data_polya_paths[overlap_data_polya_paths$Method == "mpmap_nopaths",]$Method <- "vg mpmap"

overlap_data_polya_paths$Graph <- factor(overlap_data_polya_paths$Graph, levels = c("With transcript paths", "Without transcript paths"))
plotOverlapBenchmark(overlap_data_polya_paths, wes_cols, "plots/polya_rna/vg_sim_overlap_polya_paths")


########

# 
# overlap_data_mir <- overlap_data %>%
#   filter(Reads == "ENCSR958UOC_rep1_uni") 
# 
# overlap_data_mir$Method <- recode_factor(overlap_data_mir$Method, 
#                                      "bowtie2" = "Bowtie2",
#                                      "bowtie2_vs_end" = "Bowtie2 (vs end)",
#                                      "bowtie2_vs_local" = "Bowtie2 (vs local)",
#                                      "hisat2_nosplice" = "HISAT2", 
#                                      "star_nosplice" = "STAR", 
#                                      "star_encode" = "STAR (encode)",
#                                      "map" = "vg map (def)", 
#                                      "map_fast" = "vg map", 
#                                      "mpmap" = "vg mpmap (multi)", 
#                                      "mpmap_nomulti" = "vg mpmap")
# 
# overlap_data_mir$Reads <- recode_factor(overlap_data_mir$Reads, 
#                                          "ERR187607_uni" = "Training set (old)", 
#                                          "ENCSR958UOC_rep1_uni" = "Training set (rep1)")
# 
# overlap_data_mir$Graph = recode_factor(overlap_data_mir$Graph, 
#                                         "linear" = "Linear",
#                                         "gencode100" = "Linear",
#                                         "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")
# 
# overlap_data_mir$FacetCol <- recode_factor(overlap_data_mir$Threshold, 
#                                              "Overlap >= 90%" = "Overlap >= 90%", 
#                                              "Overlap >= 10%" = "Overlap >= 10%")
# 
# overlap_data_mir$FacetRow <- recode_factor(overlap_data_mir$Filter, 
#                                              "Unfiltered" = "Unfiltered", 
#                                              "Low quality bases filtered (< 4)" = "Low quality bases filtered (< 4)")
# 
# 
# overlap_data_mir_main <- overlap_data_mir %>%
#   filter(Method != "Bowtie2 (vs end)") %>%
#   filter(Method != "Bowtie2 (vs local)") %>%
#   filter(Method != "STAR") %>% 
#   filter(Method != "vg map (def)") %>%
#   filter(Method != "vg mpmap (multi)")
#   
# wes_cols_mir_main <- wes_cols[c(6, seq(1, 5))]
# 
# plotOverlapBenchmark(overlap_data_mir_main, wes_cols_mir_main, "plots/micro_rna/vg_sim_overlap_mir_main")
# 
# 
# overlap_data_mir_train <- overlap_data_mir %>%
#   filter(Method != "vg mpmap (multi)")
# 
# wes_cols_mir_train <- c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2"))
# 
# plotOverlapBenchmark(overlap_data_mir_train, wes_cols_mir_train, "plots/micro_rna/vg_sim_overlap_mir_train")
# 
# 
