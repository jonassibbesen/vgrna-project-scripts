
# plot_vg_sim_distance.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_names = F)
  colnames(data) <- c("Count", "Distance", "IsMapped", "MapQ", "Method")

  data <- data %>%
    select(-Method) %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])

  if (grepl("dist_gamp_", basename(filename))) {

    data <- data %>%
      mutate(Method = paste(Method, "gamp", sep = "_"))
  }
    
  return(data)
}

distance_threshold <- 100

distance_data_raw <- map_dfr(list.files(path = "./methods", pattern=".*_dist_gam.*.txt", full.names = T, recursive = T), parse_file)

distance_data <- distance_data_raw %>%
  mutate(Correct = Distance <= distance_threshold) 


########


distance_data_polya <- distance_data %>%
  filter(Type == "polya_rna")

distance_data_polya$Method <- recode_factor(distance_data_polya$Method, 
                                      "hisat2" = "HISAT2",
                                      "star" = "STAR",
                                      "map" = "vg map (def)", 
                                      "map_fast" = "vg map", 
                                      "mpmap" = "vg mpmap (gam)", 
                                      "mpmap_gamp" = "vg mpmap")

distance_data_polya <- distance_data_polya %>%
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (gam)")

distance_data_polya$FacetCol <- "Simulated reads"
distance_data_polya$FacetRow <- ""


distance_data_polya_main <- distance_data_polya %>%
  filter(Graph != "1kg_nonCEU_af001_gencode100_genes")

distance_data_polya_main$Graph = recode_factor(distance_data_polya_main$Graph, 
                                        "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                        "gencode100" = "Spliced reference")

for (reads in unique(distance_data_polya_main$Reads)) {
  
  distance_data_polya_main_reads <- distance_data_polya_main %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQ(distance_data_polya_main_reads, wes_cols, paste("plots/polya_rna/vg_sim_distance_polya_main_d", distance_threshold, "_", reads, sep = ""))
}


distance_data_polya_gene <- distance_data_polya %>%
  filter(Graph != "gencode100") %>%
  filter(Method != "STAR") %>%
  filter(Method != "HISAT2")

distance_data_polya_gene$Graph = recode_factor(distance_data_polya_gene$Graph, 
                                         "1kg_nonCEU_af001_gencode100" = "Whole genome graph",
                                         "1kg_nonCEU_af001_gencode100_genes" = "Exons only graph")

for (reads in unique(distance_data_polya_gene$Reads)) {
  
  distance_data_polya_gene_reads <- distance_data_polya_gene %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQ(distance_data_polya_gene_reads, wes_cols[c(3,4)], paste("plots/polya_rna/vg_sim_distance_polya_gene_d", distance_threshold, "_", reads, sep = ""))
}


# distance_data_polya_paths <- distance_data_polya %>%
#   filter(Graph == "1kg_nonCEU_af001_gencode100" | Method == "STAR") 
# 
# distance_data_polya_paths[distance_data_polya_paths$Method == "HISAT2",]$Graph <- "Without transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "STAR",]$Graph <- "Without transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "vg map",]$Graph <- "With transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "vg mpmap (gam)",]$Graph <- "With transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "vg mpmap (gamp)",]$Graph <- "With transcript paths"
# 
# distance_data_polya_paths[distance_data_polya_paths$Method == "map_fast_nopaths",]$Graph <- "Without transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "map_fast_nopaths",]$Method <- "vg map"
# 
# distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths",]$Graph <- "Without transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths",]$Method <- "vg mpmap (gam)"
# 
# distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths_gamp",]$Graph <- "Without transcript paths"
# distance_data_polya_paths[distance_data_polya_paths$Method == "mpmap_nopaths_gamp",]$Method <- "vg mpmap (gamp)"
# 
# distance_data_polya_paths$Graph <- factor(distance_data_polya_paths$Graph, levels = c("With transcript paths", "Without transcript paths"))
# plotDistanceBenchmark(distance_data_polya_paths, wes_cols, "plots/polya_rna/vg_sim_distance_polya_paths")


########
