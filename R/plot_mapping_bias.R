
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

min_site_count = 1000
min_site_mapq = 20
frac_site_mapq = 0.075

coverage_data <- map_dfr(list.files(pattern=".*_allele_cov.txt", full.names = T, recursive = T), parse_file)

coverage_data <- coverage_data %>%
  filter(Graph != "gencode85") %>%
  filter(Graph != "1kg_nonCEU_af001_gencode85") %>%
  filter(Method != "map_fast_nopaths") %>%
  filter(Method != "mpmap_nopaths")
  
# low_mapq_sites <- coverage_data %>%
#   filter(Method != "Bowtie2 (vs end)") %>%
#   filter(Method != "Bowtie2 (vs local)") %>%
#   filter(Method != "STAR") %>% 
#   filter(Method != "vg map") %>%
#   filter(Method != "vg mpmap (multi)") %>%
#   group_by(VariantPosition, Reads, Method, Graph) %>%
#   mutate(read_count = Count * (UpReadCount + DownReadCount) / 2) %>%
#   mutate(low_mapq_read_count = ifelse(MapQ < min_site_mapq, read_count, 0)) %>%
#   summarise(count = sum(read_count), frac_low = sum(low_mapq_read_count) / sum(read_count)) %>%
#   filter(count >= min_site_count) %>%
#   filter(frac_low > frac_site_mapq) %>%
#   ungroup() %>%
#   select(VariantPosition, Reads) %>%
#   add_column(low_mapq_site = T)

coverage_data_mq <- coverage_data %>%
  #left_join(low_mapq_sites, by = c("VariantPosition", "Reads")) %>%
  # filter(is.na(low_mapq_site)) %>%
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
  mutate(len = ifelse(len > 15, 15, len)) %>%
  mutate(len = ifelse(len < -15, -15, len)) %>%
  group_by(Reads, Method, Graph, var, len) %>%
  summarise(n = n(), ref_count = sum(ref), alt_count = sum(alt), frac_mean = mean(frac))   

set.seed(123456)

pdf(paste("plots/polya_rna/vg_sim_mapping_bias_polya_main.pdf", sep = ""), height = 3, width = 7, pointsize = 12)
coverage_data_mq_polya_filt %>% 
    ggplot(aes(y = frac_mean, x = len, color = Method, linetype = Graph, shape = Graph, label = sprintf("%0.3f", round(frac_mean, digits = 3)))) +
    geom_line(size = 0.75) + 
    geom_point(data = subset(coverage_data_mq_polya_filt, len == 0), size = 1.5) +  
    geom_text_repel(data = subset(coverage_data_mq_polya_filt, len == 0), size = 2, fontface = 2, box.padding = 0.75) +  
    geom_hline(yintercept = 0.5, size = 0.25, linetype = 1, alpha = 0.75) + 
    facet_grid(rows = vars(Graph)) + 
    scale_color_manual(values = wes_cols) +
    ylim(c(0.3, 0.7)) +
    xlab("Allele length") +
    ylab("Mean fraction mapped reads to alt allele") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size = 10))
dev.off()

coverage_data_mq_polya_filt_simple <- coverage_data_mq_polya_filt %>% 
  filter(Method == "STAR" | (Method == "vg mpmap" & Graph == "1000g (no-CEU)")) 

pdf(paste("plots/polya_rna/vg_sim_mapping_bias_polya_simple.pdf", sep = ""), height = 3, width = 7, pointsize = 12)
coverage_data_mq_polya_filt_simple %>% 
  filter(Method == "STAR" | (Method == "vg mpmap" & Graph == "1000g (no-CEU)")) %>%
  ggplot(aes(y = frac_mean, x = len, color = Method, label = sprintf("%0.3f", round(frac_mean, digits = 3)))) +
  geom_line(size = 0.75) + 
  geom_point(data = subset(coverage_data_mq_polya_filt_simple, len == 0), size = 1.5) +  
  geom_text_repel(data = subset(coverage_data_mq_polya_filt_simple, len == 0), size = 2, fontface = 2, box.padding = 0.75) +  
  geom_hline(yintercept = 0.5, size = 0.25, linetype = 1, alpha = 0.75) + 
  scale_color_manual(values = wes_cols) +
  ylim(c(0.4, 0.6)) +
  xlab("Allele length") +
  ylab("Mean fraction mapped reads to alt allele") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 10))
dev.off()


coverage_data_mq_polya_binom <- coverage_data_mq_polya %>%
  mutate(ref = (ref_up + ref_down) / 2) %>%
  mutate(alt = (alt_up + alt_down) / 2) %>%
  filter(ref + alt >= min_count) %>%
  mutate(len = ifelse(len > 10, 10, len)) %>%
  mutate(len = ifelse(len < -10, -10, len)) %>%
  rowwise() %>%
  mutate(binom_test = binom.test(x = c(as.integer(ref), as.integer(alt)), alternative = "two.sided")$p.value) %>%
  group_by(Reads, Method, Graph, var, len) %>%
  summarise(n = n(), n_binom = sum(binom_test < 0.05))   

pdf(paste("plots/polya_rna/vg_sim_mapping_bias_polya_binom.pdf", sep = ""), height = 3, width = 7, pointsize = 12)
coverage_data_mq_polya_binom %>% 
  filter(Method == "STAR" | (Method == "vg mpmap" & Graph == "1000g (no-CEU)")) %>%
  ggplot(aes(y = n_binom / n, x = len, color = Method)) +
  geom_line(size = 0.75) + 
  scale_color_manual(values = wes_cols) +
  xlab("Allele length") +
  ylab("Fraction variants with binomial test p-value < 0.05") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 8))
dev.off()



########


min_count = 1000

coverage_data_mq_mir <- coverage_data_mq %>%
  filter(Reads == "ENCSR958UOC_rep1_uni")

coverage_data_mq_mir$Method <- recode_factor(coverage_data_mq_mir$Method, 
                                         "bowtie2" = "Bowtie2",
                                         "bowtie2_vs_end" = "Bowtie2 (vs end)",
                                         "bowtie2_vs_local" = "Bowtie2 (vs local)",
                                         "hisat2_nosplice" = "HISAT2", 
                                         "star_nosplice" = "STAR", 
                                         "star_encode" = "STAR (encode)",
                                         "map" = "vg map (def)", 
                                         "map_fast" = "vg map", 
                                         "mpmap" = "vg mpmap (multi)", 
                                         "mpmap_nomulti" = "vg mpmap")

coverage_data_mq_mir$Reads <- recode_factor(coverage_data_mq_mir$Reads, 
                                        "ERR187607_uni" = "Training set (old)", 
                                        "ENCSR958UOC_rep1_uni" = "Training set (rep1)")

coverage_data_mq_mir$Graph = recode_factor(coverage_data_mq_mir$Graph, 
                                        "linear" = "Linear",
                                        "gencode100" = "Linear",
                                        "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")

coverage_data_mq_mir$Reads <- NULL

coverage_data_mq_mir_filt <- coverage_data_mq_mir %>%
  mutate(ref = (ref_up + ref_down) / 2) %>%
  mutate(alt = (alt_up + alt_down) / 2) %>%
  filter(ref + alt >= min_count) %>%
  mutate(frac = alt / (ref + alt))


plotMappingBiasMir <- function(mapping_bias_data, wes_cols_mir, suffix) {

  set.seed(1234)
  
  mapping_bias_data <- mapping_bias_data %>%
    filter(var == "SNV")  
  
  p <- mapping_bias_data %>% 
    ggplot(aes(y = frac, x = Method, color = Method, shape = Graph)) +
    geom_point(size = 1.25, position = position_jitterdodge(jitter.width = 2)) +
    geom_hline(yintercept = 0.5, size = 0.2, linetype = 1, alpha = 0.5) + 
    #facet_grid(cols = vars(Graph)) + 
    scale_color_manual(values = wes_cols_mir) +
    ylim(c(0,1)) +
    xlab("") +
    ylab("Mean fraction mapped reads to alt allele (SNVs)") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size = 10))
  
  pdf(paste("plots/micro_rna/vg_sim_mapping_bias_mir_", suffix, ".pdf", sep = ""), height = 3, width = 7, pointsize = 12)
  print(p)  
  dev.off()
  
  mapping_bias_data_mean <- mapping_bias_data %>%
    group_by(Method, Graph, var, len) %>%
    summarise(n = n(), sum = sum(ref + alt), frac_mean = mean(frac), frac_sd = sd(frac)) %>%
    rowwise() %>%
    mutate(frac_sd = ifelse(frac_sd > frac_mean, frac_mean, frac_sd))
  
  mapping_bias_data_mean <- mapping_bias_data_mean %>%
    ungroup() %>%
    add_row(Method = "Bowtie2", Graph = "1000g (no-CEU)", var = "SNV", len = 0, n = 0, frac_mean = NA, frac_sd = NA) %>%
    add_row(Method = "STAR (encode)", Graph = "1000g (no-CEU)", var = "SNV", len = 0, n = 0, frac_mean = NA, frac_sd = NA)
  
  mapping_bias_data_mean$Graph = recode_factor(mapping_bias_data_mean$Graph, 
                                                       "Linear" = "Linear",
                                                       "1000g (no-CEU)" = "1000g\n(no-CEU)")
  
  print(mapping_bias_data_mean)
  
  p <- mapping_bias_data_mean %>% 
    ggplot(aes(x = Graph, y = frac_mean, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
    geom_errorbar(aes(ymin = frac_mean - frac_sd, ymax = frac_mean + frac_sd), width = 0.25, size = 0.35, position=position_dodge(0.5), alpha = 0.75) +
    geom_hline(yintercept = 0.5, size = 0.2, linetype = 1, alpha = 0.5) + 
    scale_fill_manual(values = wes_cols_mir) +
    scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
    xlab("") +
    ylab("Mean fraction mapped reads to alt allele (SNVs)") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size = 10))
  
  pdf(paste("plots/micro_rna/vg_sim_mapping_bias_mir_mean_", suffix, ".pdf", sep = ""), height = 4, width = 4, pointsize = 12)
  print(p)
  dev.off()
}


coverage_data_mq_mir_filt_main <- coverage_data_mq_mir_filt %>%
  filter(Method != "Bowtie2 (vs end)") %>%
  filter(Method != "Bowtie2 (vs local)") %>%
  filter(Method != "STAR") %>% 
  filter(Method != "vg map (def)") %>%
  filter(Method != "vg mpmap (multi)")

wes_cols_mir_main <- wes_cols[c(6, seq(1, 5))]

plotMappingBiasMir(coverage_data_mq_mir_filt_main, wes_cols_mir_main, "main_filt")


coverage_data_mq_mir_filt_train <- coverage_data_mq_mir_filt %>%
  filter(Method != "vg mpmap (multi)")

wes_cols_mir_train <- c(wes_palette("Darjeeling1"), wes_palette("Darjeeling2"))

plotMappingBiasMir(coverage_data_mq_mir_filt_train, wes_cols_mir_train, "train_filt")


coverage_data_mq_mir_filt_main %>%
  filter(Method == "Bowtie2") %>%
  filter(Graph == "Linear") %>%
  select(ref, alt, frac, Graph, var) %>%
  print(n = 100)

