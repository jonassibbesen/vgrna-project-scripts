
# plot_mapping_bias.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

#source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping_r1/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  if (grepl("mpmap", filename)) {
    
    data <- read_table2(filename, col_types = "iiiciciii")
    
  } else {
    
    data <- read_table2(filename, col_types = "iiciciii") %>%
      mutate(AllelicMapQ = MapQ)
  }
  
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = paste(dir_split[6], dir_split[7], sep = "_")) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
 
  if (grepl("sim_vg", dir_split[6])) {
    
    data <- data %>%
      add_column(Simulation = "Simulated reads (vg)")
    
  } else if (grepl("sim_rsem", dir_split[6])) {
    
    data <- data %>%
      add_column(Simulation = "Simulated reads (RSEM)")
    
  } else {
    
    stopifnot(FALSE)
  }
  
  return(data)
}

min_mapq = 30

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_allele_cov.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Graph != "gencode80") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode80") 

coverage_data_amq <- coverage_data %>%
  filter(Method == "mpmap") %>%
  mutate(Method = "mpmap_amq") %>%
  mutate(MapQ = AllelicMapQ)

coverage_data <- rbind(coverage_data, coverage_data_amq)

coverage_data_mq <- coverage_data %>%
  mutate(UpReadCount = ifelse(MapQ < min_mapq, 0, UpReadCount)) %>%
  mutate(DownReadCount = ifelse(MapQ < min_mapq, 0, DownReadCount)) %>%
  group_by(VariantPosition, AlleleId, AlleleType, RelativeAlleleLength, Type, Reads, Method, Graph, Simulation) %>%
  summarise(UpReadCount = sum(UpReadCount), DownReadCount = sum(DownReadCount)) 

coverage_data_mq <- full_join(coverage_data_mq[coverage_data_mq$AlleleId == 1,], coverage_data_mq[coverage_data_mq$AlleleId == 2,], by = c("VariantPosition", "Type", "Reads", "Method", "Graph", "Simulation"))

coverage_data_mq <- coverage_data_mq %>% 
  filter((AlleleType.x != AlleleType.y) & (AlleleType.x == "REF" | AlleleType.y == "REF")) %>% 
  mutate(ref_up = ifelse(AlleleType.x == "REF", UpReadCount.x, UpReadCount.y)) %>%
  mutate(ref_down = ifelse(AlleleType.x == "REF", DownReadCount.x, DownReadCount.y)) %>%
  mutate(alt_up = ifelse(AlleleType.x != "REF", UpReadCount.x, UpReadCount.y)) %>%
  mutate(alt_down = ifelse(AlleleType.x != "REF", DownReadCount.x, DownReadCount.y)) %>%
  mutate(var = ifelse(AlleleType.x == "REF", AlleleType.y, AlleleType.x)) %>%
  mutate(len = ifelse(AlleleType.x == "REF", RelativeAlleleLength.y, RelativeAlleleLength.x)) %>%
  filter(var != "COM") 

coverage_data_mq$var <- factor(coverage_data_mq$var, levels = c("SNV", "INS", "DEL"))

coverage_data_mq$var = recode_factor(coverage_data_mq$var, 
                                     "SNV" = "SNV", 
                                     "INS" = "Insertion", 
                                     "DEL" = "Deletion")


########


coverage_data_mq_bias <- coverage_data_mq %>%
  filter(Type == "polya_rna")

coverage_data_mq_bias$Method = recode_factor(coverage_data_mq_bias$Method, 
                                        "hisat2" = "HISAT2", 
                                        "star" = "STAR",
                                        "map" = "vg map (def)", 
                                        "map_fast" = "vg map",
                                        "mpmap" = "vg mpmap", 
                                        "mpmap_amq" = "vg mpmap (aMapQ)", 
                                        "star_wasp" = "WASP (STAR)")

coverage_data_mq_bias[coverage_data_mq_bias$Method == "WASP (STAR)",]$Graph <- "1kg_NA12878_gencode100"

coverage_data_mq_bias$FacetCol <- coverage_data_mq_bias$Simulation
coverage_data_mq_bias$FacetRow <- coverage_data_mq_bias$Graph


coverage_data_mq_bias_debug <- coverage_data_mq_bias

# for (reads in unique(coverage_data_mq_bias_debug$Reads)) {
#   
#   coverage_data_mq_bias_debug_reads <- coverage_data_mq_bias_debug %>%
#     filter(Reads == reads)
#   
#   plotMappingBiasBenchmark(coverage_data_mq_bias_debug_reads, wes_cols, paste("plots/sim_bias/vg_sim_r1_mapping_bias_debug_", reads, sep = ""))
# }


coverage_data_mq_bias_main <- coverage_data_mq_bias %>%
  filter(Method != "WASP (STAR)") %>%
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (aMapQ)") %>%
  filter(Graph != "1kg_NA12878_gencode100") %>%
  filter(Graph != "1kg_NA12878_exons_gencode100")

coverage_data_mq_bias_main$Graph = recode_factor(coverage_data_mq_bias_main$Graph,
                                            "1kg_nonCEU_af001_gencode100" = "Spliced pan-\ngenome graph",
                                            "gencode100" = "Spliced\nreference")

coverage_data_mq_bias_main$FacetRow <- coverage_data_mq_bias_main$Graph

# for (reads in unique(coverage_data_mq_bias_main$Reads)) {
#   
#   coverage_data_mq_bias_main_reads <- coverage_data_mq_bias_main %>%
#     filter(Reads == reads)
#   
#   plotMappingBiasBenchmark(coverage_data_mq_bias_main_reads, wes_cols, paste("plots/sim_bias/vg_sim_r1_mapping_bias_main_", reads, sep = ""))
# }


########

coverage_data_mq_bias_binom <- coverage_data_mq_bias %>%
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (aMapQ)") 

coverage_data_mq_bias_binom$Graph = recode_factor(coverage_data_mq_bias_binom$Graph,
                                                 "1kg_nonCEU_af001_gencode100" = "Spliced pangenome graph",
                                                 "1kg_NA12878_gencode100" = "Spliced personal graph",
                                                 "1kg_NA12878_exons_gencode100" = "Spliced personal graph",
                                                 "gencode100" = "Spliced reference")

coverage_data_mq_bias_binom$FacetCol <- coverage_data_mq_bias_binom$var
coverage_data_mq_bias_binom$FacetRow <- coverage_data_mq_bias_binom$Simulation

for (reads in unique(coverage_data_mq_bias_binom$Reads)) {
  
  print(reads)
  
  coverage_data_mq_bias_binom_reads <- coverage_data_mq_bias_binom %>%
    filter(Reads == reads)
  
  plotMappingBiasBinomBenchmark(coverage_data_mq_bias_binom_reads, wes_cols, paste("plots/sim_bias/vg_sim_r1_mapping_bias_binom_", reads, sep = ""))
}

