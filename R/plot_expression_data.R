
# plot_expression_correlation2.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("ggrepel")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quantification/")

exp_data_hap_pos_all <- list()
exp_data_hap_exp_all <- list()
exp_data_stats_all <- list()

for (f in list.files(pattern = ".*SRR11534701kg.*RData", full.names = T, recursive = T)) { 
  
  load(f)
  
  exp_data_hap_pos_all[[f]] <- exp_data_hap_pos
  exp_data_hap_exp_all[[f]] <- exp_data_hap_exp
  exp_data_stats_all[[f]] <- exp_data_stats
}

wes_cols <- c(wes_palette("Rushmore1")[c(1,3,4,5)], wes_palette("Chevalier1"))

do.call(bind_rows, exp_data_hap_exp_all) %>%
  group_by(Reads, Method, Graph) %>%
  filter(tpm_est > 0) %>%
  summarise(min_tpm_est = min(tpm_est)) %>%
  print(n = 100)


######


exp_data_hap_pos_all_roc <- do.call(bind_rows, exp_data_hap_pos_all) %>%
  group_by(HaplotypePosterior, Reads, Method, Graph) %>%
  summarise(TP_pos = sum(TP_pos),
            TN_pos = sum(TN_pos),
            FP_pos = sum(FP_pos),
            FN_pos = sum(FN_pos),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Graph) %>%
  arrange(desc(HaplotypePosterior), .by_group = T) %>%
  mutate(TP_pos_cs = cumsum(TP_pos),
         TN_pos_cs = cumsum(TN_pos),
         FP_pos_cs = cumsum(FP_pos),
         FN_pos_cs = cumsum(FN_pos),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) 

exp_data_hap_pos_all_roc_sub <- exp_data_hap_pos_all_roc %>%
  group_by(Reads, Method, Graph) %>%
  filter((row_number() - 1) %% 100 == 0 | HaplotypePosterior == 0 | HaplotypePosterior == 1)

exp_data_hap_pos_all_roc_sub$Method = recode_factor(exp_data_hap_pos_all_roc_sub$Method, 
                                                   "kallisto" = "Kallisto",
                                                   "salmon" = "Salmon",
                                                   "rpvg_exact" = "rpvg")

exp_data_hap_pos_all_roc_sub$Graph = recode_factor(exp_data_hap_pos_all_roc_sub$Graph, 
                                               "1kg_NA12878_gencode100" = "Personal",
                                               "1kg_NA12878_gencode100_genes" = "Personal",
                                               "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)",
                                               "1kg_nonCEU_af001_gencode100_genes" = "1000g (no-CEU)",
                                               "1kg_all_af001_gencode100" = "1000g (all)",
                                               "1kg_all_af001_gencode100_genes" = "1000g (all)")

pdf("posterior_rocs.pdf")

exp_data_hap_pos_all_roc_sub %>%
  mutate(TPR_count = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), FPR_count = FP_pos_cs / max(FP_pos_cs + TN_pos_cs)) %>%
  ggplot(aes(y = TPR_count, x = FPR_count, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text_repel(size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols) +
  xlab("Transcript expression error (FPR)") +
  ylab("Transcript expression sensitivity (TPR)") +
  #xlim(c(0, 0.025)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos_all_roc_sub %>%
  mutate(TPR_count = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), PPV_count = TP_pos_cs / (TP_pos_cs + FP_pos_cs)) %>%
  ggplot(aes(y = TPR_count, x = PPV_count, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text_repel(size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols) +
  xlab("Transcript expression precision (PPV)") +
  ylab("Transcript expression sensitivity (TPR)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos_all_roc_sub %>%
  ggplot(aes(y = TP_pos_cs, x = FP_pos_cs, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text_repel(size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols) +
  xlab("Number incorrect expressed transcripts") +
  ylab("Number correct expressed transcripts") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos_all_roc_sub %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) %>%
  ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Method, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  geom_text_repel(size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols) +
  xlab("Fraction TPM on not expressed transcripts") +
  ylab("Fraction TPM on expressed transcripts") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

dev.off()


######


exp_data_hap_exp_all_roc <- do.call(bind_rows, exp_data_hap_exp_all) %>%
  group_by(tpm_est, Reads, Method, Graph) %>%
  summarise(TP_pos = sum(TP_pos),
            TN_pos = sum(TN_pos),
            FP_pos = sum(FP_pos),
            FN_pos = sum(FN_pos),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Graph) %>%
  arrange(desc(tpm_est), .by_group = T) %>%
  mutate(TP_pos_cs = cumsum(TP_pos),
         TN_pos_cs = cumsum(TN_pos),
         FP_pos_cs = cumsum(FP_pos),
         FN_pos_cs = cumsum(FN_pos),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) 

exp_data_hap_exp_all_roc_sub <- exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Graph) %>%
  filter((row_number() - 1) %% 100 == 0 | tpm_est == max(tpm_est) | tpm_est == min(tpm_est))

exp_data_hap_exp_all_roc_sub$Method = recode_factor(exp_data_hap_exp_all_roc_sub$Method, 
                                                    "kallisto" = "Kallisto",
                                                    "salmon" = "Salmon",
                                                    "rpvg_exact" = "rpvg")

exp_data_hap_exp_all_roc_sub$Graph = recode_factor(exp_data_hap_exp_all_roc_sub$Graph, 
                                                   "1kg_NA12878_gencode100" = "Personal",
                                                   "1kg_NA12878_gencode100_genes" = "Personal",
                                                   "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)",
                                                   "1kg_nonCEU_af001_gencode100_genes" = "1000g (no-CEU)",
                                                   "1kg_all_af001_gencode100" = "1000g (all)",
                                                   "1kg_all_af001_gencode100_genes" = "1000g (all)")

pdf("expression_rocs.pdf")

exp_data_hap_exp_all_roc_sub %>%
  mutate(TPR_count = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), FPR_count = FP_pos_cs / max(FP_pos_cs + TN_pos_cs)) %>%
  ggplot(aes(y = TPR_count, x = FPR_count, linetype = Graph, shape = Graph, color = Method)) +
  geom_line(size = 1) +
  scale_color_manual(values = wes_cols) +
  xlab("Transcript expression error (FPR)") +
  ylab("Transcript expression sensitivity (TPR)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_exp_all_roc_sub %>%
  mutate(TPR_count = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), PPV_count = TP_pos_cs / (TP_pos_cs + FP_pos_cs)) %>%
  ggplot(aes(y = TPR_count, x = PPV_count, linetype = Graph, shape = Graph, color = Method)) +
  geom_line(size = 1) +
  scale_color_manual(values = wes_cols) +
  xlab("Transcript expression precision (PPV)") +
  ylab("Transcript expression sensitivity (TPR)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_exp_all_roc_sub %>%
  ggplot(aes(y = TP_pos_cs, x = FP_pos_cs, linetype = Graph, shape = Graph, color = Method)) +
  geom_line(size = 1) +
  scale_color_manual(values = wes_cols) +
  xlab("Number incorrect expressed transcripts") +
  ylab("Number correct expressed transcripts") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_exp_all_roc_sub %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) %>%
  ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, linetype = Graph, shape = Graph, color = Method)) +
  geom_line(size = 1) +
  scale_color_manual(values = wes_cols) +
  xlab("Fraction TPM on not expressed transcripts") +
  ylab("Fraction TPM on expressed transcripts") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

dev.off()


######


exp_data_all_stats_all_bars <- do.call(bind_rows, exp_data_all_stats_all) %>%
  filter(is_hap | is.na(is_hap))

exp_data_all_stats_all_bars$Type <- "All"
exp_data_all_stats_all_bars[!is.na(exp_data_all_stats_all_bars$is_hap),]$Type <- "NA12878"

exp_data_all_stats_all_bars$Method = recode_factor(exp_data_all_stats_all_bars$Method, 
                                                    "kallisto" = "Kallisto",
                                                    "salmon" = "Salmon",
                                                    "rpvg_exact" = "rpvg")

exp_data_all_stats_all_bars$Graph = recode_factor(exp_data_all_stats_all_bars$Graph, 
                                                  "1kg_NA12878_gencode100" = "Personal\n(NA12878)",
                                                  "1kg_NA12878_gencode100_genes" = "Personal\n(NA12878)",
                                                  "1kg_nonCEU_af001_gencode100" = "1000g\n(no-CEU)",
                                                  "1kg_nonCEU_af001_gencode100_genes" = "1000g\n(no-CEU)",
                                                  "1kg_all_af001_gencode100" = "1000g\n(all)",
                                                  "1kg_all_af001_gencode100_genes" = "1000g\n(all)")

exp_data_all_stats_all_bars$Type = recode_factor(exp_data_all_stats_all_bars$Type, 
                                                  "All" = "All transcripts",
                                                  "NA12878" = "Sample transcripts")
                                            

pdf("general_stats.pdf")

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  xlab("") +
  ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = num_hap_error, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = ExpCorrect, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Fraction correct expressed") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = Pearson, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = PearsonLog, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Log Pearson correlation") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = Spearman, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_all_stats_all_bars %>%
  ggplot(aes(x = Graph, y = ARD_mean, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(cols = vars(Type)) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Mean absolute difference") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

dev.off()


###############

# 
# pdf("expression_rocs2.pdf")
# 
# exp_data_hap_exp_all_roc_sub %>%
#   filter(Graph != "Personal") %>%
#   mutate(TPR_count = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), PPV_count = TP_pos_cs / (TP_pos_cs + FP_pos_cs)) %>%
#   ggplot(aes(y = TPR_count, x = PPV_count, linetype = Graph, shape = Graph, color = Method)) +
#   geom_line(size = 1) +
#   scale_color_manual(values = wes_cols) +
#   xlim(c(0, 1)) +
#   ylim(c(0, 1)) +
#   xlab("Haplotype-specific transcript precision") +
#   ylab("Haplotype-specific transcript sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio=1) +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size = 14))
# 
# dev.off()
# 
# pdf("expression_rocs2.pdf")
# 
# exp_data_hap_exp_all_roc_sub %>%
#   filter(Graph != "Personal") %>%
#   mutate(TPR_count = TP_pos_cs / max(TP_pos_cs + FN_pos_cs), PPV_count = TP_pos_cs / (TP_pos_cs + FP_pos_cs)) %>%
#   ggplot(aes(y = TPR_count, x = PPV_count, linetype = Graph, shape = Graph, color = Method)) +
#   geom_line(size = 1.5) +
#   scale_color_manual(values = wes_cols) +
#   xlim(c(0, 1)) +
#   ylim(c(0, 1)) +
#   xlab("Haplotype-specific transcript precision") +
#   ylab("Haplotype-specific transcript sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio=1) +
#   theme(strip.background = element_blank()) +
#   theme(legend.position="bottom", legend.box="vertical") +
#   theme(text = element_text(size = 16))
# 
# dev.off()
# 
# 
# pdf("hap_error.pdf", width = 5)
# 
# exp_data_all_stats_all_bars %>%
#   filter(Graph != "Personal\n(NA12878)") %>%
#   filter(Type == "All transcripts") %>%
#   ggplot(aes(x = Graph, y = hap_error, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
#   scale_fill_manual(values = wes_cols) +
#   xlab("") +
#   ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size = 16))
# 
# dev.off()
# 
# pdf("hap_error2.pdf", width = 4)
# 
# exp_data_all_stats_all_bars %>%
#   filter(Graph != "Personal\n(NA12878)") %>%
#   filter(Type == "All transcripts") %>%
#   ggplot(aes(x = Graph, y = hap_error, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
#   scale_fill_manual(values = wes_cols) +
#   xlab("") +
#   ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(legend.position="bottom") +
#   theme(text = element_text(size = 16))
# 
# dev.off()
# 
# 
# pdf("spearman.pdf", width = 9)
# 
# exp_data_all_stats_all_bars %>%
#   ggplot(aes(x = Graph, y = Spearman, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
#   facet_grid(cols = vars(Type)) +
#   scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
#   scale_fill_manual(values = wes_cols) +
#   xlab("") +
#   ylab("Spearman expression correlation") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size = 16))
# 
# dev.off()
# 
# pdf("spearman2.pdf", width = 7)
# 
# exp_data_all_stats_all_bars %>%
#   ggplot(aes(x = Graph, y = Spearman, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
#   facet_grid(cols = vars(Type)) +
#   scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
#   scale_fill_manual(values = wes_cols) +
#   xlab("") +
#   ylab("Spearman expression correlation") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(legend.position="bottom") +
#   theme(text = element_text(size = 15))
# 
# dev.off()
# 
# 
