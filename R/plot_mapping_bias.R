
# plot_mapping_bias.R

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
  
  data <- read_table2(filename, col_types = "iicci")
  data <- data %>%
    add_column(Reads = dir_split[6]) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8])
  
  return(data)
}

coverage_data_h1 <- map_dfr(list.files(pattern=".*_allele_cov_h1.txt", full.names = T, recursive = T), parse_file)
coverage_data_h2 <- map_dfr(list.files(pattern=".*_allele_cov_h2.txt", full.names = T, recursive = T), parse_file)

coverage_data <- full_join(coverage_data_h1, coverage_data_h2, by = c("MapQ", "AllelePosition", "Reads", "Method", "Graph"))

coverage_data <- coverage_data %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85")

coverage_data <- coverage_data %>% 
  filter((AlleleType.x != AlleleType.y) & (AlleleType.x == "REF" | AlleleType.y == "REF")) %>% 
  rowwise() %>% 
  mutate(ref = ifelse(AlleleType.x == "REF", Count.x, Count.y)) %>%
  mutate(alt = ifelse(AlleleType.x != "REF", Count.x, Count.y)) %>%
  mutate(var = ifelse(AlleleType.x == "REF", AlleleType.y, AlleleType.x)) %>%
  mutate(var_len = ifelse(AlleleType.x == "REF", RelativeAlleleLength.y, RelativeAlleleLength.x)) %>%
  filter(var != "COM") 

coverage_data$var <- factor(coverage_data$var, levels = c("SNV", "INS", "DEL"))
coverage_data$var = recode_factor(coverage_data$var, "SNV" = "SNV", "INS" = "Insertion", "DEL" = "Deletion")
  
coverage_data$Method = recode_factor(coverage_data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")
coverage_data$Graph = recode_factor(coverage_data$Graph, "gencode100" = "Spliced reference", "1kg_NA12878_exons_gencode100" = "Personal (NA12878)", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

coverage_data_mq30 <- coverage_data %>%
  filter(MapQ > 30) %>%
  #mutate(ref = ifelse(var == "Deletion", ref * 101 / (var_len + 101), ref)) %>%
  filter(ref + alt >= 10) %>%
  mutate(frac = alt / (ref + alt)) %>%
  mutate(var_len = ifelse(var_len > 20, 20, var_len)) %>%
  mutate(var_len = ifelse(var == "Deletion", -1 * var_len, var_len)) %>%
  group_by(Reads, Method, Graph, var, var_len) %>%
  summarise(ref_count = sum(ref), alt_count = sum(alt), frac_median = median(frac))   

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

pdf("rsem_sim_benchmark_bias_mp30_c10_median.pdf", height = 5, width = 9)
coverage_data_mq30 %>% 
  ggplot(aes(y = frac_median, x = var_len, color = Method)) +
  geom_line(size = 0.75) + 
  facet_grid(rows = vars(Graph)) + 
  scale_color_manual(values = wes_cols) +
  ylim(c(0.2,0.55)) +
  xlab("Allele length") +
  ylab("Median fraction mapped reads to alternative allele") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14)) 
dev.off()
