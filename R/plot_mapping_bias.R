
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
  
  data <- read_table2(filename, col_types = "iiciciii")
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  return(data)
}

coverage_data <- map_dfr(list.files(pattern=".*_allele_cov.txt", full.names = T, recursive = T), parse_file)

coverage_data_mq <- coverage_data %>%
  filter(MapQ >= 30) %>%
  group_by(VariantPosition, AlleleId, AlleleType, RelativeAlleleLength, Reads, Method, Graph) %>%
  summarise(UpReadCount = sum(UpReadCount), DownReadCount = sum(DownReadCount)) 

coverage_data_mq <- full_join(coverage_data_mq[coverage_data_mq$AlleleId == 1,], coverage_data_mq[coverage_data_mq$AlleleId == 2,], by = c("VariantPosition", "Reads", "Method", "Graph"))

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
coverage_data_mq$var = recode_factor(coverage_data_mq$var, "SNV" = "SNV", "INS" = "Insertion", "DEL" = "Deletion")
  
coverage_data_mq$Method = recode_factor(coverage_data_mq$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap")
coverage_data_mq$Graph = recode_factor(coverage_data_mq$Graph, "gencode100" = "Spliced reference", "1kg_NA12878_exons_gencode100" = "Personal (NA12878)", "1kg_NA12878_gencode100" = "Personal (NA12878)", "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

coverage_data_mq_filt <- coverage_data_mq %>%
  mutate(ref = (ref_up + ref_down) / 2) %>%
  mutate(alt = (alt_up + alt_down) / 2) %>%
  filter(ref + alt >= 20) %>%
  mutate(frac = alt / (ref + alt)) %>%
  mutate(len = ifelse(len > 20, 20, len)) %>%
  mutate(len = ifelse(len < -20, -20, len)) %>%
  group_by(Reads, Method, Graph, var, len) %>%
  summarise(ref_count = sum(ref), alt_count = sum(alt), frac_mean = mean(frac))   

coverage_data_mq_filt <- coverage_data_mq_filt %>% 
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af01_gencode100") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") 

wes_cols <- c(wes_palette("Darjeeling1")[c(1,2,3,5)])

for (reads in unique(coverage_data_mq_filt$Reads)) {
  
  png(paste("rsem_sim_benchmark_bias_", reads, "_mp30_c20_mean_methods.png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  p <- coverage_data_mq_filt %>% 
    filter(Reads == reads) %>%
    ggplot(aes(y = frac_mean, x = len, color = Method)) +
    geom_line(size = 0.75) + 
    facet_grid(rows = vars(Graph)) + 
    scale_color_manual(values = wes_cols) +
    ylim(c(0.2, 0.55)) +
    xlab("Allele length") +
    ylab("Mean fraction mapped reads to alternative allele") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14))
  print(p)
  dev.off()
  
  png(paste("rsem_sim_benchmark_bias_", reads, "_mp30_c20_mean_graphs.png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  p <- coverage_data_mq_filt %>% 
    filter(Reads == reads) %>%
    ggplot(aes(y = frac_mean, x = len, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 0.75) + 
    facet_grid(rows = vars(Method)) + 
    scale_color_manual(values = wes_cols) +
    ylim(c(0.2, 0.55)) +
    xlab("Allele length") +
    ylab("Mean fraction mapped reads to alternative allele") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14))
  print(p)
  dev.off()
}
