
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

set.seed(1234)

exp_data_hap_pos_all <- list()
exp_data_hap_exp_all <- list()
exp_data_stats_all <- list()

prepareData <- function(data) {
  
  data <- data %>%
    filter(Method != "rsem_k1k") %>%
    filter(Method != "salmon_w10k")

  data$Reads = recode_factor(data$Reads,
                             "sim_vg_SRR1153470" = "Simulated data",
                             "real_SRR1153470" = "Real data")
  
  data$Method = recode_factor(data$Method,
                              "kallisto" = "Kallisto",
                              "kallisto_strand" = "Kallisto",
                              "salmon" = "Salmon",
                              "rsem" = "RSEM",
                              "rpvg4_exact_gam" = "rpvg (gam)",
                              "rpvg4_exact" = "rpvg (gamp)",
                              "rpvg7" = "rpvg7_gibbs",
                              "rpvg7_exact" = "rpvg7",
                              "rpvg7_exact_w64" = "rpvg7_w64",
                              "rpvg7_exact_w1000" = "rpvg7_w1000",
                              "rpvg8_exact" = "rpvg8")
  
  data$Graph = recode_factor(data$Graph,
                             "1kg_NA12878_gencode100" = "NA12878",
                             "1kg_NA12878_gencode100_decoy" = "NA12878",
                             "1kg_NA12878_gencode100_genes" = "NA12878",
                             "1kg_nonCEU_af001_gencode100" = "no-CEU",
                             "1kg_nonCEU_af001_gencode100_decoy" = "no-CEU",
                             "1kg_nonCEU_af001_gencode100_genes" = "no-CEU",
                             "1kg_all_af001_gencode100" = "All",
                             "1kg_all_af001_gencode100_decoy" = "All",
                             "1kg_all_af001_gencode100_genes" = "All")
  
  data$Reads <- factor(data$Reads, levels = c("Simulated data", "Real data"))
  data$Method <- factor(data$Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg (gam)", "rpvg (gamp)", "rpvg7", "rpvg7_gibbs", "rpvg7_w64", "rpvg7_w1000", "rpvg8"))
  data$Graph <- factor(data$Graph, levels = c("NA12878", "no-CEU", "All"))
  
  return(data)
}

for (f in list.files(pattern = ".*SRR11534701kg.*RData", full.names = T, recursive = T)) { 

  if (grepl("rpvg_debug", f)) {
    
    next 
  }
  
  if (grepl("rpvg5", f)) {
    
    next 
  }
  
  if (grepl("rpvg6", f)) {
    
    next 
  }
  
  if (grepl("debug", f)) {
    
    next 
  }
  
  print(f)
  
  load(f)
  
  exp_data_hap_pos_all[[f]] <- prepareData(exp_data_hap_pos)
  exp_data_hap_exp_all[[f]] <- prepareData(exp_data_hap_exp)
  exp_data_stats_all[[f]] <- prepareData(exp_data_stats)
}

do.call(bind_rows, exp_data_hap_exp_all) %>%
  group_by(Reads, Method, Graph) %>%
  filter(tpm_est > 0) %>%
  summarise(min_tpm_est = min(tpm_est), min_tpm_est = min(tpm_est)) %>%
  print(n = 100)


######


exp_data_hap_pos_all_roc <- do.call(bind_rows, exp_data_hap_pos_all) %>%
  group_by(HaplotypePosterior, Reads, Method, Graph) %>%
  summarise(TP = sum(TP),
            TN = sum(TN),
            FP = sum(FP),
            FN = sum(FN),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Graph) %>%
  arrange(desc(HaplotypePosterior), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         TN_cs = cumsum(TN),
         FP_cs = cumsum(FP),
         FN_cs = cumsum(FN),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) %>%
  mutate(TPR_count = TP_cs / max(TP_cs + FN_cs), FPR_count = FP_cs / max(FP_cs + TN_cs)) %>%
  mutate(PPV_count = TP_cs / (TP_cs + FP_cs)) %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs))
  
exp_data_hap_pos_all_roc_points <- exp_data_hap_pos_all_roc %>% 
  mutate(HaplotypePosterior_floor = floor(HaplotypePosterior * 5) / 5) %>%
  group_by(HaplotypePosterior_floor, Reads, Method, Graph) %>%
  mutate(is_min = (HaplotypePosterior == min(HaplotypePosterior))) %>%
  filter(is_min) %>%
  mutate(HaplotypePosterior = HaplotypePosterior_floor)

pdf("plots/posterior_rocs_debug.pdf", width = 9, height = 9)

exp_data_hap_pos_all_roc %>%
  ggplot(aes(y = TPR_count, x = FPR_count, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_pos_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_pos_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  xlab("Transcript expression error (FPR)") +
  ylab("Transcript expression sensitivity (TPR)") +
  #xlim(c(0, 0.025)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos_all_roc %>%
  ggplot(aes(y = TPR_count, x = PPV_count, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_pos_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_pos_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  xlab("Transcript expression precision (PPV)") +
  ylab("Transcript expression sensitivity (TPR)") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos_all_roc %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_pos_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_pos_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  xlab("Number incorrect expressed transcripts") +
  ylab("Number correct expressed transcripts") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_pos_all_roc %>%
  ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Method, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_pos_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_pos_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  xlab("Fraction TPM on not expressed transcripts") +
  ylab("Fraction TPM on expressed transcripts") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

dev.off()


######


exp_data_hap_exp_all_roc <- do.call(bind_rows, exp_data_hap_exp_all) %>%
  mutate(tpm_est_sig = signif(tpm_est, 1)) %>%
  mutate(tpm_est_sig = ifelse(tpm_est_sig >= 1, 1, tpm_est_sig)) %>%
  group_by(tpm_est_sig, Reads, Method, Graph) %>%
  summarise(TP = sum(TP),
            TN = sum(TN),
            FP = sum(FP),
            FN = sum(FN),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Graph) %>%
  arrange(desc(tpm_est_sig), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         TN_cs = cumsum(TN),
         FP_cs = cumsum(FP),
         FN_cs = cumsum(FN),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) %>%
  mutate(TPR_count = TP_cs / max(TP_cs + FN_cs), FPR_count = FP_cs / max(FP_cs + TN_cs)) %>%
  mutate(PPV_count = TP_cs / (TP_cs + FP_cs)) %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) %>%
  filter(tpm_est_sig > 0)

exp_data_hap_exp_all_roc_points <- exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Graph) %>%
  mutate(is_min = (tpm_est_sig == min(tpm_est_sig))) %>%
  filter(is_min | tpm_est_sig == 1 | tpm_est_sig == 0.01 | tpm_est_sig == 0.0001) 


bla <- exp_data_hap_exp_all_roc %>% filter(Graph == "no-CEU") %>% filter(Method == "Kallisto") %>% filter(Reads == "real_SRR1153470")


pdf("plots/expression_rocs_debug.pdf")

exp_data_hap_exp_all_roc %>%
  ggplot(aes(y = TPR_count, x = FPR_count, linetype = Graph, shape = Graph, color = Method, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_exp_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  xlab("Transcript expression error (FPR)") +
  ylab("Transcript expression sensitivity (TPR)") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_exp_all_roc %>%
  ggplot(aes(y = TPR_count, x = PPV_count, linetype = Graph, shape = Graph, color = Method, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_exp_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  xlab("Transcript expression precision (PPV)") +
  ylab("Transcript expression sensitivity (TPR)") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_exp_all_roc %>%
  ggplot(aes(y = TP_cs, x = FP_cs, linetype = Graph, shape = Graph, color = Method, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_exp_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  xlab("Number incorrect expressed transcripts") +
  ylab("Number correct expressed transcripts") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

exp_data_hap_exp_all_roc %>%
  ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, linetype = Graph, shape = Graph, color = Method, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points, size = 1.5) +
  geom_text_repel(data = exp_data_hap_exp_all_roc_points, size = 2, fontface = 2) +
  facet_grid(rows = vars(Reads)) +
  xlab("Fraction TPM on not expressed transcripts") +
  ylab("Fraction TPM on expressed transcripts") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))

dev.off()


######


exp_data_stats_all_bars <- do.call(bind_rows, exp_data_stats_all) %>%
  filter(Type != "Transcript")

pdf("plots/general_stats_debug.pdf")

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Graph, y = frac_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  xlab("") +
  ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Graph, y = frac_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  xlab("") +
  ylab("Fraction estimated count for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Graph, y = num_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Graph, y = num_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (count)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Graph, y = ExpCorrect, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Fraction correct expressed") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  ggplot(aes(x = Graph, y = Pearson_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  ggplot(aes(x = Graph, y = Pearson_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation (count)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  ggplot(aes(x = Graph, y = Spearman_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  ggplot(aes(x = Graph, y = Spearman_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation (count)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  ggplot(aes(x = Graph, y = ARD_mean_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Mean absolute difference (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  ggplot(aes(x = Graph, y = ARD_mean_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type + Truncated, scales="free") +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Mean absolute difference (count)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

dev.off()


###############

wes_cols <- wes_palette("Darjeeling2")

exp_data_hap_pos_all_roc_points_sim <- exp_data_hap_pos_all_roc_points %>%
  filter(Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Simulated data")

pdf("plots/sim_posterior_roc_pre.pdf", pointsize = 12)

exp_data_hap_pos_all_roc %>%
  filter(Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Simulated data") %>%
  ggplot(aes(y = TPR_count, x = PPV_count, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_pos_all_roc_points_sim, size = 2) +
  geom_text_repel(data = exp_data_hap_pos_all_roc_points_sim, size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols[c(4,5)]) +
  scale_shape_discrete(name = "Transcripts") +
  scale_linetype_discrete(name = "Transcripts") +
  facet_grid(cols = vars(Reads)) +
  coord_fixed() +
  xlab("Transcript expression precision") +
  ylab("Transcript expression sensitivity") +
  guides(linetype = FALSE, shape = FALSE) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 14))

dev.off()


exp_data_hap_pos_all_roc_points_real <- exp_data_hap_pos_all_roc_points %>%
  filter(Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Real data") %>%
  filter(Graph != "NA12878") 

pdf("plots/real_posterior_roc_pre.pdf", pointsize = 12)

exp_data_hap_pos_all_roc %>%
  filter(Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Real data") %>%
  ungroup() %>%
  add_row(HaplotypePosterior = 0, Reads = "Real data", Method = "rpvg (gamp)", Graph = "NA12878") %>%
  mutate(Graph = factor(Graph, levels = c("NA12878", "no-CEU", "All"))) %>%
  mutate(Method = factor(Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg (gam)", "rpvg (gamp)"))) %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Graph, shape = Graph, label = HaplotypePosterior)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_pos_all_roc_points_real, size = 2) +
  geom_text_repel(data = exp_data_hap_pos_all_roc_points_real, size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols[c(4,5)]) +
  scale_shape_discrete(name = "Transcripts") +
  scale_linetype_discrete(name = "Transcripts") +
  facet_grid(cols = vars(Reads)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  xlab("Number of transcripts expressed on non-NA12878 haplotypes") +
  ylab("Number of transcripts expressed on NA12878 haplotypes") +
  annotation_logticks() +
  guides(color = FALSE) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 14))

dev.off()


exp_data_hap_exp_all_roc_points_sim <- exp_data_hap_exp_all_roc_points %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Simulated data")

pdf("plots/sim_expression_roc_pre.pdf", pointsize = 12)

exp_data_hap_exp_all_roc %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Simulated data") %>%
  ggplot(aes(y = TPR_count, x = PPV_count, color = Method, linetype = Graph, shape = Graph, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_sim, size = 2) +
  geom_text_repel(data = exp_data_hap_exp_all_roc_points_sim, size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols) +
  scale_shape_discrete(name = "Transcripts") +
  scale_linetype_discrete(name = "Transcripts") +
  facet_grid(cols = vars(Reads)) +
  coord_fixed() +
  xlab("Transcript expression precision") +
  ylab("Transcript expression sensitivity") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 14))

dev.off()


exp_data_hap_exp_all_roc_points_real <- exp_data_hap_exp_all_roc_points %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Graph != "NA12878") %>%
  filter(Reads == "Real data")

pdf("plots/real_expression_roc_pre.pdf", pointsize = 12)

exp_data_hap_exp_all_roc %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg (gam)" | Method == "rpvg (gamp)") %>%
  filter(Reads == "Real data") %>%
  ungroup() %>%
  add_row(tpm_est_sig = 0, Reads = "Real data", Method = "rpvg (gamp)", Graph = "NA12878") %>%
  mutate(Graph = factor(Graph, levels = c("NA12878", "no-CEU", "All"))) %>%
  mutate(Method = factor(Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg (gam)", "rpvg (gamp)"))) %>%
  ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Graph, shape = Graph, label = tpm_est_sig)) +
  geom_line(size = 1) +
  geom_point(data = exp_data_hap_exp_all_roc_points_real, size = 2) +
  geom_text_repel(data = exp_data_hap_exp_all_roc_points_real, size = 3, fontface = 2) +
  scale_color_manual(values = wes_cols[c(1,2,4,5)]) +
  scale_shape_discrete(name = "Transcripts") +
  scale_linetype_discrete(name = "Transcripts") +
  facet_grid(cols = vars(Reads)) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("Number of transcripts expressed on non-NA12878 haplotypes") +
  ylab("Number of transcripts expressed on NA12878 haplotypes") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 14))

dev.off()


exp_data_stats_all_bars <- exp_data_stats_all_bars %>%
  filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg (gam)" | Method == "rpvg (gamp)")

exp_data_stats_all_bars$Graph = recode_factor(exp_data_stats_all_bars$Graph,
                                              "NA12878" = "NA12878",
                                              "no-CEU" = "no-CEU",
                                              "All" = "All")

pdf("plots/sim_real_hap_tpm_error_tpm.pdf", width = 6, pointsize = 12)

exp_data_stats_all_bars %>%
  filter(Graph != "NA12878") %>%
  filter(Type == "All") %>%
  filter(!Truncated) %>%
  ggplot(aes(x = Graph, y = frac_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  facet_grid(cols = vars(Reads)) +
  scale_fill_manual(values = wes_cols[c(1,2,4,5)]) +
  xlab("") +
  ylab("Fraction TPM on non-NA12878 haplotypes") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16))

dev.off()


pdf("plots/sim_real_hap_tpm_error_count.pdf", width = 6, pointsize = 12)

exp_data_stats_all_bars %>%
  filter(Graph != "NA12878") %>%
  filter(Type == "All") %>%
  filter(!Truncated) %>%
  ggplot(aes(x = Graph, y = frac_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  facet_grid(cols = vars(Reads)) +
  scale_fill_manual(values = wes_cols[c(1,2,4,5)]) +
  xlab("") +
  ylab("Fraction counts on non-NA12878 haplotypes") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16))

dev.off()


exp_data_stats_all_bars_sim <- exp_data_stats_all_bars %>%
  filter(Reads == "Simulated data") %>%
  filter(Type != "Transcript")  %>%
  filter(Truncated)

exp_data_stats_all_bars_sim <- exp_data_stats_all_bars_sim %>%
  ungroup() %>%
  add_row(Reads = "Simulated data", Method = "RSEM", Graph = "no-CEU", Type = "All", Spearman_tpm = 0, Spearman_count = 0, ARD_mean_tpm = 0, ARD_mean_count = 0) %>%
  add_row(Reads = "Simulated data", Method = "RSEM", Graph = "no-CEU", Type = "Haplotype", Spearman_tpm = 0, Spearman_count = 0, ARD_mean_tpm = 0, ARD_mean_count = 0) %>%
  add_row(Reads = "Simulated data", Method = "RSEM", Graph = "All", Type = "All", Spearman_tpm = 0, Spearman_count = 0, ARD_mean_tpm = 0, ARD_mean_count = 0) %>%
  add_row(Reads = "Simulated data", Method = "RSEM", Graph = "All", Type = "Haplotype", Spearman_tpm = 0, Spearman_count = 0, ARD_mean_tpm = 0, ARD_mean_count = 0)

exp_data_stats_all_bars_sim$Type <- as.factor(exp_data_stats_all_bars_sim$Type)
exp_data_stats_all_bars_sim$Method <- factor(exp_data_stats_all_bars_sim$Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg (gam)", "rpvg (gamp)"))
exp_data_stats_all_bars_sim$Graph <- factor(exp_data_stats_all_bars_sim$Graph, levels = c("NA12878", "no-CEU", "All"))

pdf("plots/sim_spearman_all_hap_tpm.pdf", width = 6, pointsize = 12)

ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "All"), aes(x = Graph, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "Haplotype"), aes(x = Graph, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.5), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(cols = vars(Reads)) +
  xlab("") +
  ylab("Expression Spearman correlation (TPM)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16))

dev.off()

pdf("plots/sim_spearman_all_hap_count.pdf", width = 6, pointsize = 12)

ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "All"), aes(x = Graph, y = Spearman_count, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "Haplotype"), aes(x = Graph, y = Spearman_count, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.5), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(cols = vars(Reads)) +
  xlab("") +
  ylab("Expression Spearman correlation (counts)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16))

dev.off()


pdf("plots/sim_mard_all_hap_tpm.pdf", width = 6, pointsize = 12)

ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "All"), aes(x = Graph, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "Haplotype"), aes(x = Graph, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.5), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(cols = vars(Reads)) +
  xlab("") +
  ylab("Mean absolute relative expression difference (TPM)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16))

dev.off()

pdf("plots/sim_mard_all_hap_count.pdf", width = 6, pointsize = 12)

ggplot() +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "All"), aes(x = Graph, y = ARD_mean_count, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = subset(exp_data_stats_all_bars_sim, Type == "Haplotype"), aes(x = Graph, y = ARD_mean_count, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
  scale_alpha_manual(name = "Transcripts", values = c(1, 0.5), labels = c("All", "NA12878"), drop = F) +
  scale_fill_manual(values = wes_cols) +
  facet_grid(cols = vars(Reads)) +
  xlab("") +
  ylab("Mean absolute relative expression difference (counts)") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size = 16))

dev.off()


exp_data_stats_all_bars %>%
  filter(Graph != "NA12878") %>%
  filter(Type == "All") %>%
  filter(Truncated) %>%
  select(Reads, Method, Graph, n, frac_hap_error_count, frac_hap_error_tpm) %>%
  print(n = 100)

exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Graph) %>%
  filter(tpm_est_sig == 0.01) %>%
  summarise(max_FP_cs = max(FP_cs)) %>%
  print(n = 100)
