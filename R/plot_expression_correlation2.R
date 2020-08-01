
# plot_expression_correlation2.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quantification/")

parse_exp_data_hap_pos <- function(filename) {
  
  load(filename)
  data <- exp_data_hap_pos
  
  return(data)
}

parse_exp_data_ic_count <- function(filename) {
  
  load(filename)
  data <- exp_data_ic_count
  
  return(data)
}

parse_exp_data_all_stats <- function(filename) {
  
  load(filename)
  data <- exp_data_all_stats
  
  return(data)
}

exp_data_hap_pos <- map_dfr(list.files(pattern = ".*_af001_gencode100.*.RData", full.names = T, recursive = T), parse_exp_data_hap_pos)
exp_data_ic_count <- map_dfr(list.files(pattern = ".*_af001_gencode100.*.RData", full.names = T, recursive = T), parse_exp_data_ic_count)
exp_data_all_stats <- map_dfr(list.files(pattern = ".*_af001_gencode100.*.RData", full.names = T, recursive = T), parse_exp_data_all_stats)

exp_data_hap_pos$Graph = recode_factor(exp_data_hap_pos$Graph, "1kg_all_af001_gencode100_genes" = "1000g (genes)", "1kg_all_af001_gencode100" = "1000g", "1kg_nonCEU_af001_gencode100_genes" = "1000g (genes)", "1kg_nonCEU_af001_gencode100" = "1000g")
exp_data_ic_count$Graph = recode_factor(exp_data_ic_count$Graph, "1kg_all_af001_gencode100_genes" = "1000g (genes)", "1kg_all_af001_gencode100" = "1000g", "1kg_nonCEU_af001_gencode100_genes" = "1000g (genes)", "1kg_nonCEU_af001_gencode100" = "1000g")
exp_data_all_stats$Graph = recode_factor(exp_data_all_stats$Graph, "1kg_all_af001_gencode100_genes" = "1000g (genes)", "1kg_all_af001_gencode100" = "1000g", "1kg_nonCEU_af001_gencode100_genes" = "1000g (genes)", "1kg_nonCEU_af001_gencode100" = "1000g")

wes_cols <- c(wes_palette("Rushmore1"), wes_palette("Chevalier1"))

exp_data_hap_pos <- exp_data_hap_pos %>%
  group_by(Type, Reads, Method, Graph) %>%
  arrange(desc(HaplotypePosterior), .by_group = T) %>%
  mutate(TP_hap_cs = cumsum(TP_hap), FP_hap_cs = cumsum(FP_hap), 
         TP_hap_tpm_cs = cumsum(TP_hap_tpm), FP_hap_tpm_cs = cumsum(FP_hap_tpm), 
         TP_sim_cs = cumsum(TP_sim), FP_sim_cs = cumsum(FP_sim), 
         TP_sim_tpm_cs = cumsum(TP_sim_tpm), FP_sim_tpm_cs = cumsum(FP_sim_tpm),
         TP_pos_cs = cumsum(TP_pos),
         TN_pos_cs = cumsum(TN_pos),
         FP_pos_cs = cumsum(FP_pos),
         FN_pos_cs = cumsum(FN_pos)) %>%
  filter((row_number() - 1) %% 100 == 0)
  

pdf("hap.pdf")

exp_data_hap_pos %>%
  mutate(Sensitivity = TP_hap_cs / max(TP_hap_cs), Error = FP_hap_cs / max(TP_hap_cs + FP_hap_cs)) %>%
  ggplot(aes(y = Sensitivity, x = Error, color = Method, linetype = Graph, shape = Graph)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text(aes(label = HaplotypePosterior), hjust = 1, vjust = 0) +
  facet_grid(row = vars(Type)) +
  scale_color_manual(values = wes_cols) +
  #xlim(c(0, 0.05)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos %>%
  mutate(Sensitivity = TP_hap_tpm_cs / max(TP_hap_tpm_cs), Error = FP_hap_tpm_cs / max(TP_hap_tpm_cs + FP_hap_tpm_cs)) %>%
  ggplot(aes(y = Sensitivity, x = Error, color = Method, linetype = Graph, shape = Graph)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text(aes(label = HaplotypePosterior), hjust = 1, vjust = 0) + 
  facet_grid(row = vars(Type)) +
  scale_color_manual(values = wes_cols) +
  #coord_fixed() +
  #xlim(c(0, 0.1)) +
  #ylim(c(0.9, 1)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos %>%
  mutate(Sensitivity = TP_sim_cs / max(TP_sim_cs), Error = FP_sim_cs / max(TP_sim_cs + FP_sim_cs)) %>%
  ggplot(aes(y = Sensitivity, x = Error, color = Method, linetype = Graph, shape = Graph)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text(aes(label = HaplotypePosterior), hjust = 1, vjust = 0) + 
  facet_grid(row = vars(Type)) +
  scale_color_manual(values = wes_cols) +
 # xlim(c(0, 0.05)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos %>%
  mutate(Sensitivity = TP_sim_tpm_cs / max(TP_sim_tpm_cs), Error = FP_sim_tpm_cs / max(TP_sim_tpm_cs + FP_sim_tpm_cs)) %>%
  ggplot(aes(y = Sensitivity, x = Error, color = Method, linetype = Graph, shape = Graph)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text(aes(label = HaplotypePosterior), hjust = 1, vjust = 0) + 
  facet_grid(row = vars(Type)) +
  scale_color_manual(values = wes_cols) +
  #coord_fixed() +
  #xlim(c(0, 0.1)) +
  #ylim(c(0.9, 1)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos %>%
  mutate(TPR = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), FPR = FP_pos_cs / max(FP_pos_cs + TN_pos_cs)) %>%
  ggplot(aes(y = TPR, x = FPR, color = Method, linetype = Graph, shape = Graph)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text(aes(label = HaplotypePosterior), hjust = 1, vjust = 0) + 
  facet_grid(row = vars(Type)) +
  scale_color_manual(values = wes_cols) +
  #xlim(c(0, 0.025)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

dev.off()

#exp_data$Method = recode_factor(exp_data$Method, "kallisto" = "Kallisto", "salmon" = "Salmon", "rpvg" = "rpvg")
#exp_data$Graph = recode_factor(exp_data$Graph, "1kg_NA12878_gencode100" = "Personal\n(NA12878)", "1kg_NA12878_gencode100_genes" = "Personal\n(NA12878)", "1kg_all_af001_gencode100" = "1000g\n(all)", "1kg_all_af001_gencode100_genes" = "1000g\n(all)", "1kg_nonCEU_af001_gencode100" = "1000g\n(no-CEU)", "1kg_nonCEU_af001_gencode100_genes" = "1000g\n(no-CEU)")

exp_data_all_stats <- exp_data_all_stats %>%
  replace_na(list(is_hap = "ALL")) %>%
  filter(is_hap != F)

pdf("exp_rsem_sim_hap_error_frac.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_hap_error_num.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = num_hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

exp_data_ic_count <- exp_data_ic_count %>%
  filter(!is_hap & tpm_est > 0) %>%
  add_column(count = 1) %>%
  group_by(tpm_est, Type, Reads, Method, Graph) %>%
  summarise(count = sum(count))

pdf("exp_rsem_sim_hap_error_cs.pdf")
exp_data_ic_count %>% 
  group_by(Type, Reads, Method, Graph) %>%
  arrange(desc(tpm_est), .by_group = T) %>%
  mutate(count_cs = cumsum(count)) %>%
  ggplot(aes(y = count_cs, x = tpm_est, color = Method, linetype = Graph)) +
  geom_line(size = 1.5) + 
  geom_vline(xintercept = 0.0001, linetype = "dashed") +
  facet_grid(cols = vars(Type)) +
  scale_x_continuous(trans='log10') +
  scale_color_manual(values = wes_cols) +
  xlab("Estimated TPM threshold") +
  ylab("Number of transcripts on incorrect\nhaplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_hap_error_cs_log.pdf")
exp_data_ic_count %>% 
  group_by(Type, Reads, Method, Graph) %>%
  arrange(desc(tpm_est), .by_group = T) %>%
  mutate(count_cs = cumsum(count)) %>%
  ggplot(aes(y = count_cs, x = tpm_est, color = Method, linetype = Graph)) +
  geom_line(size = 1.5) + 
  geom_vline(xintercept = 0.0001, linetype = "dashed") +
  facet_grid(cols = vars(Type)) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = wes_cols) +
  xlab("Estimated TPM threshold") +
  ylab("Number of transcripts on incorrect\nhaplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_pearson.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = Pearson, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(is_hap ~ Type) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation between estimated\nand simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_pearson_log.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = PearsonLog, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(is_hap ~ Type) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation between estimated\nand simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_spearman.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = Spearman, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(is_hap ~ Type) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation between estimated\nand simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_mean_ard.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = 1 - ARD_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(is_hap ~ Type, scales = "free") +
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("1 - MARD between estimated and simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

# pdf("exp_rsem_sim_mean_ard_violin.pdf", width = 8)
# exp_data_all %>%
#   ggplot(aes(x = Graph, y = ard, fill = Method)) +
#   geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
#   scale_fill_manual(values = wes_cols) +
#   xlab("") +
#   ylab("MARD between estimated and simulated TPM") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size=18))
# dev.off()
# 
# pdf("exp_rsem_sim_mean_ard_violin_trunc.pdf", width = 8)
# exp_data_trunc %>%
#   ggplot(aes(x = Graph, y = ard, fill = Method)) +
#   geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
#   scale_fill_manual(values = wes_cols) +
#   xlab("") +
#   ylab("MARD between estimated and simulated TPM") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size=18))
# dev.off()

