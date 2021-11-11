
# plot_vg_sim_overlap.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

#source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r1/")

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


overlap_data_raw_h1 <- map_dfr(list.files(path = "./methods", pattern=".*_ovl3_vg.*h1.txt.gz", full.names = T, recursive = T), parse_file) %>%
  select(-SubstitutionBP2, -IndelBP2) %>%
  rename(SubstitutionBP = SubstitutionBP1) %>%
  rename(IndelBP = IndelBP1)

overlap_data_raw_h2 <- map_dfr(list.files(path = "./methods", pattern=".*_ovl3_vg.*h2.txt.gz", full.names = T, recursive = T), parse_file)  %>%
  select(-SubstitutionBP1, -IndelBP1) %>%
  rename(SubstitutionBP = SubstitutionBP2) %>%
  rename(IndelBP = IndelBP2)

overlap_data <- rbind(overlap_data_raw_h1, overlap_data_raw_h2)  %>%
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

#overlap_data_polya <- overlap_data_polya %>%
#  filter(Method != "vg map (def)")

overlap_data_polya$FacetCol <- "Simulated reads"
overlap_data_polya$FacetRow <- "All reads"

overlap_data_polya <- overlap_data_polya %>%
  filter(Filter == "Low quality bases filtered")


overlap_data_polya_main <- overlap_data_polya %>%
  filter(Graph != "gencode80") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode80")

overlap_data_polya_main$Graph = recode_factor(overlap_data_polya_main$Graph, 
                                         "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                         "1kg_NA12878_gencode100" = "Personal reference graph",
                                         "1kg_NA12878_exons_gencode100" = "Personal reference graph",
                                         "gencode100" = "Spliced reference")

for (reads in unique(overlap_data_polya_main$Reads)) {
  
  overlap_data_polya_main_reads <- overlap_data_polya_main %>%
    filter(Reads == reads)

  plotOverlapBenchmarkMapQ(overlap_data_polya_main_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_main_ovl", overlap_threshold, "_", reads, sep = ""))
}


for (reads in unique(overlap_data_polya_main$Reads)) {
  
  overlap_data_polya_main_reads <- overlap_data_polya_main %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQVar(overlap_data_polya_main_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_var_ovl", overlap_threshold, "_", reads, sep = ""))
}


overlap_data_polya_main_snv1 <- overlap_data_polya_main %>%
  filter(SubstitutionBP >= 1 & SubstitutionBP <= 3 & IndelBP == 0) %>%
  mutate(FacetRow = "1-3 SNVs (no indels)")

for (reads in unique(overlap_data_polya_main_snv1$Reads)) {
  
  overlap_data_polya_main_snv1_reads <- overlap_data_polya_main_snv1 %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQ(overlap_data_polya_main_snv1_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_main_snv1_ovl", overlap_threshold, "_", reads, sep = ""))
}


overlap_data_polya_main_snv4 <- overlap_data_polya_main %>%
  filter(SubstitutionBP > 3 & IndelBP == 0) %>%
  mutate(FacetRow = ">3 SNVs (no indels)")

for (reads in unique(overlap_data_polya_main_snv4$Reads)) {
  
  overlap_data_polya_main_snv4_reads <- overlap_data_polya_main_snv4 %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQ(overlap_data_polya_main_snv4_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_main_snv4_ovl", overlap_threshold, "_", reads, sep = ""))
}


overlap_data_polya_main_indel <- overlap_data_polya_main %>%
  filter(IndelBP > 0) %>%
  mutate(FacetRow = ">0 indels")

for (reads in unique(overlap_data_polya_main_indel$Reads)) {
  
  overlap_data_polya_main_indel_reads <- overlap_data_polya_main_indel %>%
    filter(Reads == reads)
  
  plotOverlapBenchmarkMapQ(overlap_data_polya_main_indel_reads, wes_cols, paste("plots/polya_rna/vg_sim_overlap_polya_main_indel_ovl", overlap_threshold, "_", reads, sep = ""))
}


overlap_data_polya_sj <- overlap_data_polya %>%
  filter(Graph != "1kg_NA12878_gencode100") %>%
  filter(Graph != "1kg_NA12878_exons_gencode100") %>%
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




########
