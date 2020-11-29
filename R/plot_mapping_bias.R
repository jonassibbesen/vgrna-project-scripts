
# plot_mapping_bias.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

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

min_mapq = 30

coverage_data <- map_dfr(list.files(path = "./methods", pattern=".*_allele_cov.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Graph != "gencode80") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode80") 

coverage_data_mq <- coverage_data %>%
  mutate(UpReadCount = ifelse(MapQ < min_mapq, 0, UpReadCount)) %>%
  mutate(DownReadCount = ifelse(MapQ < min_mapq, 0, DownReadCount)) %>%
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

coverage_data_mq$var = recode_factor(coverage_data_mq$var, 
                                     "SNV" = "SNV", 
                                     "INS" = "Insertion", 
                                     "DEL" = "Deletion")


########

min_count = 20

coverage_data_mq_polya <- coverage_data_mq %>%
  filter(Reads == data_set1) 

coverage_data_mq_polya$Method = recode_factor(coverage_data_mq_polya$Method, 
                                        "hisat2" = "HISAT2", 
                                        "star" = "STAR", 
                                        "map" = "vg map (def)", 
                                        "map_fast" = "vg map",
                                        "mpmap" = "vg mpmap")

coverage_data_mq_polya <- coverage_data_mq_polya %>%
  filter(Method != "vg map (def)")

coverage_data_mq_polya$Graph = recode_factor(coverage_data_mq_polya$Graph, 
                                        "gencode100" = "Spliced reference",
                                        "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

coverage_data_mq_polya_filt <- coverage_data_mq_polya %>%
  mutate(ref = (ref_up + ref_down) / 2) %>%
  mutate(alt = (alt_up + alt_down) / 2) %>%
  filter(ref + alt >= min_count) %>%
  mutate(frac = alt / (ref + alt)) %>%
  mutate(len = ifelse(len > 15, 16, len)) %>%
  mutate(len = ifelse(len < -15, -16, len)) %>%
  group_by(Reads, Method, Graph, var, len) %>%
  summarise(n = n(), ref_count = sum(ref), alt_count = sum(alt), frac_mean = mean(frac))   

set.seed(123)

pdf(paste("plots/polya_rna/vg_sim_mapping_bias_polya_main.pdf", sep = ""), height = 3.5, width = 6, pointsize = 12)
coverage_data_mq_polya_filt %>% 
    ggplot(aes(y = frac_mean, x = len, color = Method, linetype = Graph, shape = Graph, label = sprintf("%0.3f", round(frac_mean, digits = 3)))) +
    geom_line(size = 1) + 
    geom_point(data = subset(coverage_data_mq_polya_filt, len == 0), size = 2) +  
    geom_text_repel(data = subset(coverage_data_mq_polya_filt, len == 0), size = 3, fontface = 2, box.padding = 0.75) +  
    geom_hline(yintercept = 0.5, size = 0.5, linetype = 1, alpha = 0.75) + 
    facet_grid(rows = vars(Graph)) + 
    scale_color_manual(values = wes_cols) +
    scale_x_continuous(breaks=c(-16, -10, -5, 0, 5, 10, 16), labels = c("<-15", "-10", "-5", "SNV", "5", "10", ">15")) +
    ylim(c(0.3, 0.7)) +
    xlab("Allele length") +
    ylab("Mean fraction alt allele reads") +
    guides(color = FALSE) +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size = 13))
dev.off()

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
