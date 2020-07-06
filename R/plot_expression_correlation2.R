
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

parse_exp_data <- function(filename) {
  
  load(filename)
  data <- exp_data_sim
  
  return(data)
}

exp_data <- map_dfr(list.files(pattern = ".*.RData", full.names = T, recursive = T), parse_exp_data)

exp_data$Method = recode_factor(exp_data$Method, "kallisto" = "Kallisto", "salmon" = "Salmon", "rpvg" = "rpvg")
exp_data$Graph = recode_factor(exp_data$Graph, "1kg_NA12878_gencode100" = "Personal\n(NA12878)", "1kg_NA12878_gencode100_genes" = "Personal\n(NA12878)", "1kg_all_af001_gencode100" = "1000g\n(all)", "1kg_all_af001_gencode100_genes" = "1000g\n(all)", "1kg_nonCEU_af001_gencode100" = "1000g\n(no-CEU)", "1kg_nonCEU_af001_gencode100_genes" = "1000g\n(no-CEU)")

exp_data_all_stats <- exp_data %>%
  mutate(ard = abs(tpm_est - tpm_sim) / (tpm_est + tpm_sim)) %>%
  replace_na(list(ard = 0)) %>%
  mutate(non_hap_tpm = tpm_est * !is_hap) %>%
  group_by(Reads, Method, Graph) %>%
  summarise(n = n(), Pearson = cor(tpm_est, tpm_sim, method = "pearson"), Spearman = cor(tpm_est, tpm_sim, method = "spearman"), ARD_mean = mean(ard), ARD_median = median(ard), hap_error = sum(non_hap_tpm) / sum(tpm_est), num_hap_error = sum(!is_hap & tpm_est > 0))

exp_data_trunc_stats <- exp_data %>%
  mutate(tpm_est = ifelse(tpm_est < 0.0001, 0, tpm_est)) %>%
  mutate(tpm_sim = ifelse(tpm_sim < 0.0001, 0, tpm_sim)) %>%
  mutate(ard = abs(tpm_est - tpm_sim) / (tpm_est + tpm_sim)) %>%
  replace_na(list(ard = 0)) %>%
  mutate(non_hap_tpm = tpm_est * !is_hap)  %>%
  group_by(Reads, Method, Graph) %>%
  summarise(n = n(), Pearson = cor(tpm_est, tpm_sim, method = "pearson"), Spearman = cor(tpm_est, tpm_sim, method = "spearman"), ARD_mean = mean(ard), ARD_median = median(ard), hap_error = sum(non_hap_tpm) / sum(tpm_est), num_hap_error = sum(!is_hap & tpm_est > 0))

wes_cols <- c(wes_palette("Rushmore1")[c(1,3,4)])

pdf("exp_rsem_sim_hap_error_frac.pdf", width = 6)
exp_data_all_stats %>%
  ggplot(aes(x = Graph, y = hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 0.25), oob = rescale_none) +
  xlab("") +
  ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_hap_error_frac_trunc.pdf", width = 6)
exp_data_trunc_stats %>%
  ggplot(aes(x = Graph, y = hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 0.25), oob = rescale_none) +
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
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_hap_error_num_trunc.pdf", width = 6)
exp_data_trunc_stats %>%
  ggplot(aes(x = Graph, y = num_hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

exp_data2 <- exp_data %>%
  filter(!is_hap & tpm_est > 0) %>%
  filter(Graph == "1000g\n(all)") %>%
  add_column(count = 1) 

exp_data2$Graph = recode_factor(exp_data2$Graph, "1000g\n(all)" = "1000g (all)")

pdf("exp_rsem_sim_hap_error_cs.pdf")
exp_data2 %>% 
  group_by(Method, Graph) %>%
  arrange(desc(tpm_est), .by_group = T) %>%
  mutate(count_cs = cumsum(count)) %>%
  ggplot(aes(y = count_cs, x = tpm_est, color = Method, linetype = Graph)) +
  geom_line(size = 1.5) + 
  geom_vline(xintercept = 0.0001, linetype = "dashed") +
  scale_x_continuous(trans='log10') +
  scale_color_manual(values = wes_cols) +
  xlab("Estimated TPM threshold") +
  ylab("Number of transcripts on incorrect\nhaplotypes (non-NA12878)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_hap_error_cs_log.pdf")
exp_data2 %>% 
  group_by(Method, Graph) %>%
  arrange(desc(tpm_est), .by_group = T) %>%
  mutate(count_cs = cumsum(count)) %>%
  ggplot(aes(y = count_cs, x = tpm_est, color = Method, linetype = Graph)) +
  geom_line(size = 1.5) + 
  geom_vline(xintercept = 0.0001, linetype = "dashed") +
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
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation between estimated\nand simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_pearson_trunc.pdf", width = 6)
exp_data_trunc_stats %>%
  ggplot(aes(x = Graph, y = Pearson, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
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
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation between estimated\nand simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_spearman_trunc.pdf", width = 6)
exp_data_trunc_stats %>%
  ggplot(aes(x = Graph, y = Spearman, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
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
  ggplot(aes(x = Graph, y = ARD_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 0.25), oob = rescale_none) +
  xlab("") +
  ylab("MARD between estimated and simulated TPM") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=18))
dev.off()

pdf("exp_rsem_sim_mean_ard_trunc.pdf", width = 6)
exp_data_trunc_stats %>%
  ggplot(aes(x = Graph, y = ARD_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  scale_fill_manual(values = wes_cols) +
  scale_y_continuous(limits=c(0, 0.25), oob = rescale_none) +
  xlab("") +
  ylab("MARD between estimated and simulated TPM") +
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

