
# plot_vg_sim_distance.R

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
  
  data <- read_table2(filename, col_names = F)
  colnames(data) <- c("Count", "Distance", "IsMapped", "MapQ", "Method")

  data <- data %>%
    select(-Method) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])

  if (grepl("dist_gamp_", basename(filename))) {

    data <- data %>%
      mutate(Method = paste(Method, "gamp", sep = "_"))
  }
    
  return(data)
}

distance_data_raw <- map_dfr(list.files(pattern=".*_dist_gam.*.txt", full.names = T, recursive = T), parse_file)

distance_data_d2 <- distance_data_raw %>%
  mutate(Correct = Distance <= 2) %>%
  add_column(Threshold = "Distance <= 2")

distance_data_d10 <- distance_data_raw %>%
  mutate(Correct = Distance <= 10) %>%
  add_column(Threshold = "Distance <= 10")

distance_data_d20 <- distance_data_raw %>%
  mutate(Correct = Distance <= 20) %>%
  add_column(Threshold = "Distance <= 20")

distance_data_d90 <- distance_data_raw %>%
  mutate(Correct = Distance <= 90) %>%
  add_column(Threshold = "Distance <= 90")

distance_data <- bind_rows(distance_data_d2, distance_data_d10, distance_data_d20, distance_data_d90)


########


distance_data_polya <- distance_data %>%
  filter(Reads == data_set1) %>%
  filter(Threshold == "Distance <= 10" | Threshold == "Distance <= 90") 

distance_data_polya$Method <- recode_factor(distance_data_polya$Method, 
                                      "hisat2" = "HISAT2",
                                      "star" = "STAR",
                                      "map" = "vg map (def)", 
                                      "map_fast" = "vg map", 
                                      "mpmap" = "vg mpmap (gam)", 
                                      "mpmap_gamp" = "vg mpmap (gamp)")

distance_data_polya <- distance_data_polya %>%
  filter(Method != "vg map (def)")
  
distance_data_polya$Reads <- recode_factor(distance_data_polya$Reads, 
                                    "SRR1153470_uni" = "Training set", 
                                    "ENCSR000AED_rep1_uni" = "Test set")

distance_data_polya$FacetCol <- distance_data_polya$Threshold
distance_data_polya$FacetRow <- ""


distance_data_polya_main <- distance_data_polya %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") %>%
  filter(Method != "map_fast_nopaths") %>%
  filter(Method != "mpmap_nopaths") %>%
  filter(Method != "mpmap_nopaths_gamp") 

distance_data_polya_main$Graph = recode_factor(distance_data_polya_main$Graph, 
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

plotDistanceBenchmark(distance_data_polya_main, wes_cols, "plots/polya_rna/vg_sim_distance_polya_main")


distance_data_polya_sj <- distance_data_polya %>%
  filter(Graph != "gencode100" | Method == "STAR") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes") %>%
  filter(Method != "map_fast_nopaths") %>%
  filter(Method != "mpmap_nopaths") %>%
  filter(Method != "mpmap_nopaths_gamp") 

distance_data_polya_sj$Graph = recode_factor(distance_data_polya_sj$Graph, 
                                       "gencode100" = "All transcripts", 
                                       "1kg_nonCEU_af001_gencode100" = "All transcripts", 
                                       "gencode85" = "85% transcripts", 
                                       "1kg_nonCEU_af001_gencode85" = "85% transcripts")

plotDistanceBenchmark(distance_data_polya_sj, wes_cols, "plots/polya_rna/vg_sim_distance_polya_sj")


distance_data_polya_gene <- distance_data_polya %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "gencode100" | Method == "STAR") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Method != "map_fast_nopaths") %>%
  filter(Method != "mpmap_nopaths") %>%
  filter(Method != "mpmap_nopaths_gamp") 

distance_data_polya_gene$Graph = recode_factor(distance_data_polya_gene$Graph, 
                                         "gencode100" = "Whole genome", 
                                         "1kg_nonCEU_af001_gencode100" = "Whole genome",
                                         "1kg_nonCEU_af001_gencode100_genes" = "Exons only")

plotDistanceBenchmark(distance_data_polya_gene, wes_cols, "plots/polya_rna/vg_sim_distance_polya_gene")


distance_data_polya_paths <- distance_data_polya %>%
  filter(Graph == "1kg_nonCEU_af001_gencode100" | Method == "STAR") 

distance_data_polya_paths[distance_data_polya_paths$Method == "HISAT2",]$Graph <- "Without transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "STAR",]$Graph <- "Without transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "vg map",]$Graph <- "With transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "vg mpmap (gam)",]$Graph <- "With transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "vg mpmap (gamp)",]$Graph <- "With transcript paths"

distance_data_polya_paths[distance_data_polya_paths$Method == "map_fast_nopaths",]$Graph <- "Without transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "map_fast_nopaths",]$Method <- "vg map"

distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths",]$Graph <- "Without transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths",]$Method <- "vg mpmap (gam)"

distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths_gamp",]$Graph <- "Without transcript paths"
distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths_gamp",]$Method <- "vg mpmap (gamp)"

distance_data_polya_paths$Graph <- factor(distance_data_polya_paths$Graph, levels = c("With transcript paths", "Without transcript paths"))
plotDistanceBenchmark(distance_data_polya_paths, wes_cols, "plots/polya_rna/vg_sim_distance_polya_paths")


########


distance_data_mir <- distance_data %>%
  filter(Reads == "ENCSR958UOC_rep1_uni") %>%
  filter(Threshold == "Distance <= 2" | Threshold == "Distance <= 20")

distance_data_mir$Method <- recode_factor(distance_data_mir$Method, 
                                            "bowtie2" = "Bowtie2",
                                            "bowtie2_vs_end" = "Bowtie2 (vs end)",
                                            "bowtie2_vs_local" = "Bowtie2 (vs local)",
                                            "hisat2_nosplice" = "HISAT2", 
                                            "star_nosplice" = "STAR", 
                                            "star_encode" = "STAR (encode)",
                                            "map" = "vg map (def)", 
                                            "map_fast" = "vg map", 
                                            "mpmap" = "vg mpmap (multi, gam)", 
                                            "mpmap_nomulti" = "vg mpmap (gam)", 
                                            "mpmap_gamp" = "vg mpmap (gamp)",
                                            "mpmap_nomulti_gamp" = "vg mpmap (nomulti, gamp)")

distance_data_mir$Reads <- recode_factor(distance_data_mir$Reads, 
                                           "ERR187607_uni" = "Training set (old)", 
                                           "ENCSR958UOC_rep1_uni" = "Training set (rep1)")

distance_data_mir$Graph = recode_factor(distance_data_mir$Graph, 
                                             "linear" = "Linear",
                                             "gencode100" = "Linear",
                                             "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

distance_data_mir$FacetCol <- distance_data_mir$Threshold
distance_data_mir$FacetRow <- ""


distance_data_mir_main <- distance_data_mir %>%
  filter(Method != "Bowtie2 (vs end)") %>%
  filter(Method != "Bowtie2 (vs local)") %>%
  filter(Method != "STAR") %>%
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (multi, gam)") %>%
  filter(Method != "vg mpmap (nomulti, gamp)")

wes_cols_mir_main <- wes_cols[c(6, seq(1, 5))]

plotDistanceBenchmark(distance_data_mir_main, wes_cols_mir_main, "plots/micro_rna/vg_sim_distance_mir_main")

distance_data_mir_train <- distance_data_mir %>%
  filter(Method != "vg mpmap (multi, gam)") %>%
  filter(Method != "vg mpmap (nomulti, gamp)")

wes_cols_mir_train <- c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2"))

plotDistanceBenchmark(distance_data_mir_train, wes_cols_mir_train, "plots/micro_rna/vg_sim_distance_mir_train")

