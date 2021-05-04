
# parse_haploid_expression_data.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

parse_rpvg <- function(filename) {
  
  data <- read_table2(gzfile(filename), col_names = T) %>%
    group_by(ClusterID) %>%
    mutate(n_clust = n()) %>%
    filter(HaplotypeProbability >= 0.8) %>%
    filter(ReadCount > 0) %>%
    ungroup() %>%
    select(-ClusterID, -EffectiveLength) %>% 
    rename(tpm_est = TPM, count_est = ReadCount) %>%
    separate(Name, c("transcript", "hap_id"), "_") 
    
  if ("ClusterRelativeExpression" %in% names(data)) {
    
    data <- data %>%
      select(-ClusterRelativeExpression)     
  }
  
  return(data)
}

parse_salmon <- function(filename) {
  
  data <- read_table2(gzfile(filename), col_names = T) %>%
    select(-EffectiveLength) %>%
    filter(NumReads > 0) %>%
    add_column(HaplotypeProbability = 1) %>%
    add_column(n_clust = 1) %>%
    rename(tpm_est = TPM, count_est = NumReads) %>%
    separate(Name, c("transcript", "hap_id"), "_") 
  
  return(data)
}

parse_kallisto <- function(filename) {
  
  data <- read_table2(gzfile(filename), col_names = T) %>%
    select(-eff_length) %>%
    filter(est_counts > 0) %>%  
    add_column(HaplotypeProbability = 1) %>%
    add_column(n_clust = 1) %>%
    rename(Name = target_id, tpm_est = tpm, count_est = est_counts, Length = length) %>%
    separate(Name, c("transcript", "hap_id"), "_") 
  
  return(data)
}

# salmon <- parse_salmon("methods/salmon/expression/polya_rna/real/CHM13_rep1/salmon/1kg_all_af001_gencode100_decoy/salmon_1kg_all_af001_gencode100_decoy_real_CHM13_rep1/quant.sf.gz") %>%
#   add_column(Method = "salmon")
# 
# kallisto <- parse_kallisto("methods/kallisto/expression/polya_rna/real/CHM13_rep1/kallisto_strand/1kg_all_af001_gencode100/kallisto_strand_1kg_all_af001_gencode100_real_CHM13_rep1/abundance.tsv.gz") %>%
#   add_column(Method = "kallisto_strand")

rpvg <- parse_rpvg("methods/rpvg/expression/polya_rna/real/CHM13_rep1/rpvg_strand/1kg_all_af001_gencode100_unidi/rpvg_strand_mpmap_1kg_all_af001_gencode100_unidi_real_CHM13_rep1.txt.gz")  %>%
  add_column(Method = "rpvg_strand")

rpvg_v2 <- parse_rpvg("methods/rpvg/expression/polya_rna/real/CHM13_rep1/rpvg_strand/1kg_all_af001_gencode100_v2_unidi/rpvg_strand_mpmap_1kg_all_af001_gencode100_v2_unidi_real_CHM13_rep1.txt.gz")  %>%
  add_column(Method = "rpvg_strand_v2")

# hap_exp_data <- rbind(salmon, kallisto, rpvg)

hap_exp_data <- rbind(rpvg, rpvg_v2)

hap_exp_data <- hap_exp_data %>%
  group_by(Method) %>%
  mutate(count_est = count_est / sum(count_est) * 1000000) %>%
  group_by(transcript, Method) %>%
  add_column(count = 1) %>%
  arrange(desc(tpm_est), .by_group = T) %>%
  summarise(tpm_est = tpm_est, count_est = count_est, Length = Length, n_clust = n_clust, n_cs = cumsum(count)) %>% 
  mutate(major = (n_cs == 1)) 


wes_cols <- c(wes_palette("GrandBudapest1")[1], wes_palette("GrandBudapest2")[4], wes_palette("Chevalier1")[1])

# hap_exp_data$Method = recode_factor(hap_exp_data$Method,
#                             "kallisto_strand" = "Kallisto",
#                             "salmon" = "Salmon",
#                             "rpvg_strand" = "rpvg"
#                             )
# 
# hap_exp_data$Method <- factor(hap_exp_data$Method, levels = c("Kallisto", "Salmon", "rpvg"))

hap_exp_data <- hap_exp_data %>%
  add_column(Pantranscriptome = "Whole") %>%
  add_column(FacetRow = "") %>%
  add_column(FacetCol = "Real reads")

hap_exp_data_stats <- hap_exp_data %>% 
  group_by(major, Method, Pantranscriptome, FacetRow, FacetCol) %>% 
  summarise(n = n(), tpm_sum = sum(tpm_est), count_sum = sum(count_est)) %>% 
  group_by(Method, Pantranscriptome, FacetRow, FacetCol) %>% 
  mutate(n_frac = n / sum(n), tpm_sum_frac = tpm_sum / sum(tpm_sum), count_sum_frac = count_sum / sum(count_sum)) 

hap_exp_data_stats %>% print(n = 100)


pdf("real_t2t_major_count_frac_tpm.pdf", height = 4, width = 4, pointsize = 12)
hap_exp_data_stats %>%
  filter(major = TRUE) %>%
  ggplot(aes(x = Pantranscriptome, y = tpm_sum_frac, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ FacetCol) +
  xlab("") +
  ylab("Fraction TPM on major transcripts") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 13))
dev.off()

pdf("real_t2t_major_count_frac_count.pdf", height = 4, width = 4, pointsize = 12)
hap_exp_data_stats %>%
  filter(major = TRUE) %>%
  ggplot(aes(x = Pantranscriptome, y = count_sum_frac, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(FacetRow ~ FacetCol) +
  xlab("") +
  ylab("Fraction CPM on major transcripts") +
  scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 13))
dev.off()


hap_exp_data_roc_tpm <- hap_exp_data %>%
  group_by(tpm_est, Method, Pantranscriptome, FacetRow, FacetCol) %>%
  summarise(TP = sum(major == TRUE),
            FP = sum(major == FALSE)) %>%
  mutate(tpm_est_sig = signif(tpm_est, 1)) %>%
  mutate(tpm_est_sig = ifelse(tpm_est_sig >= 10, 10, tpm_est_sig)) %>%
  group_by(tpm_est_sig, Method, Pantranscriptome, FacetRow, FacetCol) %>%
  summarise(TP = sum(TP),
            FP = sum(FP)) %>%
  group_by(Method, Pantranscriptome, FacetRow, FacetCol) %>%
  arrange(desc(tpm_est_sig), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         FP_cs = cumsum(FP)) %>%
  filter(tpm_est_sig > 0)

hap_exp_data_roc_tpm_points <- hap_exp_data_roc_tpm %>%
  group_by(Method, Pantranscriptome, FacetRow, FacetCol) %>%
  mutate(is_min = (tpm_est_sig == min(tpm_est_sig))) %>%
  mutate(is_max = (tpm_est_sig == max(tpm_est_sig))) %>%
  filter(is_min | tpm_est_sig == 0.1 | tpm_est_sig == 1 | is_max) 

set.seed(1234)

pdf("real_t2t_major_minor_exp_roc_tpm.pdf", height = 5, width = 7, pointsize = 12)
hap_exp_data_roc_tpm %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = hap_exp_data_roc_tpm_points, size = 2) +
  geom_label_repel(data = hap_exp_data_roc_tpm_points[hap_exp_data_roc_tpm_points$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
  facet_grid(FacetRow ~ FacetCol) +
  scale_color_manual(values = wes_cols) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("Number of minor expressed transcripts") +
  ylab("Number of major expressed transcripts") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(legend.key.width = unit(1, "cm")) +
  theme(text = element_text(size = 15))
dev.off()




# bla %>% filter(transcript == "ENST00000346234.6") 
# bla %>% filter(transcript == "ENST00000461096.6") 
# bla %>% filter(transcript == "ENST00000317897.4") 
# bla %>% filter(transcript == "ENST00000594159.1") 
# bla %>% filter(transcript == "ENST00000396062.3") 
# bla %>% filter(transcript == "ENST00000296677.4") 
# bla %>% filter(transcript == "ENST00000568280.1") 
# bla %>% filter(transcript == "ENST00000370206.8") 
# bla %>% filter(transcript == "ENST00000379375.5") 
# bla %>% filter(transcript == "ENST00000275766.1") 
# bla %>% filter(transcript == "ENST00000317811.5") 
# bla %>% filter(transcript == "ENST00000378045.4") 
# bla %>% filter(transcript == "ENST00000329235.6") 
# bla %>% filter(transcript == "ENST00000596580.2") 
# bla %>% filter(transcript == "ENST00000368847.4") 
# 
# bla %>% filter(transcript == "ENST00000580018.3") 
# bla %>% filter(transcript == "ENST00000374259.7") 
# bla %>% filter(transcript == "ENST00000325307.11") 
# bla %>% filter(transcript == "ENST00000216252.3") 
# bla %>% filter(transcript == "ENST00000271638.2") 
# 
# bla %>% filter(transcript == "ENST00000253788.11") 
# bla %>% filter(transcript == "ENST00000323699.8") 
# bla %>% filter(transcript == "ENST00000378119.8") 
# bla %>% filter(transcript == "ENST00000228506.7") 
# bla %>% filter(transcript == "ENST00000340913.10") 
# bla %>% filter(transcript == "ENST00000592588.6") 
# bla %>% filter(transcript == "ENST00000584577.5") 
# bla %>% filter(transcript == "ENST00000223641.4") 
# bla %>% filter(transcript == "ENST00000591776.5") 
# bla %>% filter(transcript == "ENST00000221975.6") 
# 
# bla %>% filter(transcript == "ENST00000215909.9") 
# bla %>% filter(transcript == "ENST00000225964.9") 
# bla %>% filter(transcript == "ENST00000646664.1") 
# bla %>% filter(transcript == "ENST00000394077.7") 
# bla %>% filter(transcript == "ENST00000375820.9") 

