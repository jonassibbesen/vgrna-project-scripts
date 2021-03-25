
# plot_mapping_bias.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_types = "iiciciii")
  data <- data %>%
    add_column(Type = dir_split[5]) %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  return(data)
}

min_mapq = 30

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_allele_cov.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Graph != "gencode80") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode80") 

coverage_data_mq <- coverage_data %>%
  mutate(UpReadCount = ifelse(MapQ < min_mapq, 0, UpReadCount)) %>%
  mutate(DownReadCount = ifelse(MapQ < min_mapq, 0, DownReadCount)) %>%
  group_by(VariantPosition, AlleleId, AlleleType, RelativeAlleleLength, Type, Reads, Method, Graph) %>%
  summarise(UpReadCount = sum(UpReadCount), DownReadCount = sum(DownReadCount)) 

coverage_data_mq <- full_join(coverage_data_mq[coverage_data_mq$AlleleId == 1,], coverage_data_mq[coverage_data_mq$AlleleId == 2,], by = c("VariantPosition", "Type", "Reads", "Method", "Graph"))

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


min_count = 20

coverage_data_mq_polya <- coverage_data_mq %>%
  filter(Type == "polya_rna")

coverage_data_mq_polya$Method = recode_factor(coverage_data_mq_polya$Method, 
                                        "hisat2" = "HISAT2", 
                                        "star" = "STAR", 
                                        "map" = "vg map (def)", 
                                        "map_fast" = "vg map",
                                        "mpmap" = "vg mpmap")

coverage_data_mq_polya <- coverage_data_mq_polya %>%
  filter(Method != "vg map (def)")

coverage_data_mq_polya$Graph = recode_factor(coverage_data_mq_polya$Graph,
                                        "1kg_nonCEU_af001_gencode100" = "Spliced pan-\ngenome graph",
                                        "gencode100" = "Spliced reference")

coverage_data_mq_polya$FacetCol <- "Simulated reads"
coverage_data_mq_polya$FacetRow <- coverage_data_mq_polya$Graph

for (reads in unique(coverage_data_mq_polya$Reads)) {
  
  coverage_data_mq_polya_reads <- coverage_data_mq_polya %>%
    filter(Reads == reads)
  
  plotMappingBiasBenchmark(coverage_data_mq_polya_reads, wes_cols, paste("plots/polya_rna/vg_sim_mapping_bias_polya_main_", reads, sep = ""))
}


# 
# coverage_data_mq_polya_binom <- coverage_data_mq_polya %>%
#   mutate(ref = (ref_up + ref_down) / 2) %>%
#   mutate(alt = (alt_up + alt_down) / 2) %>%
#   filter(ref + alt >= min_count) %>%
#   mutate(len = ifelse(len > 10, 10, len)) %>%
#   mutate(len = ifelse(len < -10, -10, len)) %>%
#   rowwise() %>%
#   mutate(binom_test = binom.test(x = c(as.integer(ref), as.integer(alt)), alternative = "two.sided")$p.value) %>%
#   group_by(Reads, Method, Graph, var, len) %>%
#   summarise(n = n(), n_binom = sum(binom_test < 0.05))   
# 
# pdf(paste("plots/polya_rna/vg_sim_mapping_bias_polya_binom.pdf", sep = ""), height = 3, width = 7, pointsize = 12)
# coverage_data_mq_polya_binom %>% 
#   filter(Method == "STAR" | (Method == "vg mpmap" & Graph == "1000g (no-CEU)")) %>%
#   ggplot(aes(y = n_binom / n, x = len, color = Method)) +
#   geom_line(size = 0.75) + 
#   scale_color_manual(values = wes_cols) +
#   xlab("Allele length") +
#   ylab("Fraction variants with binomial test p-value < 0.05") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size = 8))
# dev.off()
# 
