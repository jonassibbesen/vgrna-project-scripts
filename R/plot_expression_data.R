
# plot_expression_data.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("ggrepel")
library("wesanderson")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quant_r1/")

set.seed(4321)

exp_data_hap_prob_all <- list()
exp_data_hap_exp_all <- list()
exp_data_stats_all <- list()

prepareData <- function(data) {
  
  data$Reads = recode_factor(data$Reads,
                              "sim_vg_SRR1153470" = "Simulated reads",
                              "real_SRR1153470" = "Real reads",
                              "sim_vg_ENCSR000AED_rep1" = "Simulated reads",
                              "real_ENCSR000AED_rep1" = "Real reads",
                              "sim_vg_r1_SRR1153470" = "Simulated reads",
                              "real_r1_SRR1153470" = "Real reads",
                              "sim_vg_r1_ENCSR000AED_rep1" = "Simulated reads",
                              "real_r1_ENCSR000AED_rep1" = "Real reads"
                             )

  data$Method = recode_factor(data$Method,
                              "rpvg" = "rpvg_r1",
                              "rpvg_strand" = "rpvg_r1",
                              "rpvg2" = "rpvg_r1_2",
                              "rpvg2_strand" = "rpvg_r1_2",
                              "rpvg3" = "rpvg_r1_3",
                              "rpvg3_strand" = "rpvg_r1_3", 
                              "rpvg4" = "rpvg_r1_4",
                              "rpvg4_strand" = "rpvg_r1_4", 
                              "rpvg5" = "rpvg_r1_5",
                              "rpvg5_strand" = "rpvg_r1_5", 
                              "rpvg6" = "rpvg_r1_6",
                              "rpvg6_strand" = "rpvg_r1_6", 
                              "rpvg7" = "rpvg_r1_7",
                              "rpvg7_strand" = "rpvg_r1_7", 
                              "rpvg8" = "rpvg_r1_8",
                              "rpvg8_strand" = "rpvg_r1_8", 
                              "rpvg9" = "rpvg_r1_9",
                              "rpvg9_strand" = "rpvg_r1_9",
                              "rpvg10" = "rpvg_r1_10",
                              "rpvg10_strand" = "rpvg_r1_10",
                              "rpvg11" = "rpvg_r1_11",
                              "rpvg11_strand" = "rpvg_r1_11",
                              "rpvg11_amq" = "rpvg_r1_11_amq",
                              "rpvg11_strand_amq" = "rpvg_r1_11_amq",
                              "rpvg12_amq" = "rpvg_r1_12_amq",
                              "rpvg12_strand_amq" = "rpvg_r1_12_amq",
                              "rpvg13_amq" = "rpvg_r1_13_amq",
                              "rpvg13_strand_amq" = "rpvg_r1_13_amq"
                             )
     
  # data$Method = recode_factor(data$Method,
  #                             "kallisto" = "Kallisto",
  #                             "kallisto_boot100" = "Kallisto (b100)",
  #                             "kallisto_strand" = "Kallisto",
  #                             "kallisto_strand_bias" = "Kallisto (bias)",
  #                             "salmon" = "Salmon",
  #                             "salmon_bias" = "Salmon (bias)",
  #                             "salmon_w5k" = "Salmon (w5000)",
  #                             "salmon_em" = "Salmon (EM)",
  #                             "salmon_vbp01" = "Salmon (VB prior 0.1)",
  #                             "salmon_vbp1" = "Salmon (VB prior 1)",
  #                             "salmon_gibbs100" = "Salmon (Gibbs 100)",
  #                             "rsem" = "RSEM (k200)",
  #                             "rsem_strand" = "RSEM (k200)",
  #                             "rsem_k1k" = "RSEM",
  #                             "rsem_k1k_pme" = "RSEM (pme)",
  #                             "rsem_strand_k1k" = "RSEM",
  #                             "rsem_k2k" = "RSEM (k2000)",
  #                             "rsem_strand_k2k" = "RSEM (k2000)",
  #                             "rpvg" = "rpvg",
  #                             "rpvg_gam" = "rpvg (single-path)",
  #                             "rpvg_strand" = "rpvg",
  #                             "rpvg_strand_gam" = "rpvg (single-path)"
  #                             )
  
  data$Graph = recode_factor(data$Graph,
                             "1kg_NA12878_gencode100" = "Sample-specific (NA12878)",
                             "1kg_NA12878_gencode100_genes" = "Sample-specific (NA12878)",
                             "1kg_EURnonCEU_af002_gencode100" = "Europe (excl. CEU)",
                             "1kg_EURnonCEU_af002_gencode100_genes" = "Europe (excl. CEU)",
                             "1kg_EURnonCEU_af002_gencode100_unidi" = "Europe (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100" = "Whole (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100_genes" ="Whole (excl. CEU)",
                             "1kg_nonCEU_af001_gencode100_unidi" = "Whole (excl. CEU)",
                             "1kg_all_af001_gencode100" = "Whole",
                             "1kg_all_af001_gencode100_genes" = "Whole",
                             "1kg_all_af001_gencode100_unidi" = "Whole",
                             "1kg_all_af001_mt_gencode100" = "Whole (mt)",
                             "1kg_all_af001_mt_gencode100_genes" = "Whole (mt)",
                             "1kg_all_af001_mt_gencode100_unidi" = "Whole (mt)",
                             "1kg_all_af001_mt_gencode100_genes_main" = "Whole (mt main)",
                             "1kg_all_af001_mt_gencode100_unidi_main" = "Whole (mt main)"
)

  data$Reads <- factor(data$Reads, levels = c("Simulated reads", "Real reads"))
#  data$Method <- factor(data$Method, levels = c("Kallisto", "Kallisto (b100)", "Kallisto (bias)", "Salmon", "Salmon (bias)", "Salmon (w5000)", "Salmon (EM)", "Salmon (VB prior 0.1)", "Salmon (VB prior 1)", "Salmon (Gibbs 100)", "RSEM", "RSEM (pme)", "RSEM (k200)", "RSEM (k2000)", "rpvg", "rpvg (single-path)"))
  data$Method <- factor(data$Method, levels = c("rpvg_paper", "rpvg_r1", "rpvg_r1_2", "rpvg_r1_3", "rpvg_r1_4", "rpvg_r1_5", "rpvg_r1_6", "rpvg_r1_7", "rpvg_r1_8", "rpvg_r1_9", "rpvg_r1_10", "rpvg_r1_11", "rpvg_r1_11_amq", "rpvg_r1_12_amq", "rpvg_r1_13_amq"))
  data$Graph <- factor(data$Graph, levels = c("Sample-specific (NA12878)", "Europe (excl. CEU)", "Whole (excl. CEU)", "Whole", "Whole (mt)", "Whole (mt main)"))
  
  data <- data %>%
    rename(Pantranscriptome = Graph)

  return(data)
}

dataset <- "SRR1153470"
#dataset <- "ENCSR000AED_rep1"

for (f in list.files(path = "./rdata", pattern = paste(".*", dataset, "1kg.*RData", sep = ""), full.names = T, recursive = T)) { 
  
  if (!grepl("rpvg", f)) {
    
    next 
  }
  
  if (grepl("rpvg2", f)) {
    
    next 
  }
  
  if (grepl("rpvg3", f)) {
    
    next 
  }
  
  if (grepl("rpvg4", f)) {
    
    next 
  }
  
  if (grepl("rpvg5", f)) {
    
    next 
  }
  
  if (grepl("rpvg6", f)) {
    
    next 
  }
  
  if (grepl("rpvg7", f)) {
    
    next 
  }
  
  if (grepl("rpvg8", f)) {
    
    next 
  }
  
  if (grepl("rpvg9", f)) {
    
    next 
  }
  
  print(f)
  load(f)
    
  exp_data_hap_prob_all[[f]] <- prepareData(exp_data_hap_prob)
  exp_data_hap_exp_all[[f]] <- prepareData(exp_data_hap_exp)
  exp_data_stats_all[[f]] <- prepareData(exp_data_stats)
}

for (f in list.files(path = "../quant/rdata", pattern = paste(".*", dataset, "1kg.*RData", sep = ""), full.names = T, recursive = T)) { 
  
  if (!grepl("rpvg", f)) {
    
    next 
  }
  
  if (grepl("1kg_nonCEU_af001_gencode100", f)) {
    
    next 
  }
  
  if (grepl("v2", f)) {
    
    next 
  }

  print(f)
  load(f)
  
  exp_data_hap_prob$Method = recode_factor(exp_data_hap_prob$Method,
                                           "rpvg" = "rpvg_paper",
                                           "rpvg_strand" = "rpvg_paper")
  
  exp_data_hap_exp$Method = recode_factor(exp_data_hap_exp$Method,
                                          "rpvg" = "rpvg_paper",
                                          "rpvg_strand" = "rpvg_paper")
  
  exp_data_stats$Method = recode_factor(exp_data_stats$Method,
                                        "rpvg" = "rpvg_paper",
                                        "rpvg_strand" = "rpvg_paper")  
  
  exp_data_hap_prob_all[[f]] <- prepareData(exp_data_hap_prob)
  exp_data_hap_exp_all[[f]] <- prepareData(exp_data_hap_exp)
  exp_data_stats_all[[f]] <- prepareData(exp_data_stats)
}

do.call(bind_rows, exp_data_hap_exp_all) %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  filter(tpm_est > 0) %>%
  summarise(min_tpm_est = min(tpm_est), max_tpm_est = max(tpm_est)) %>%
  print(n = 100)


######


exp_data_hap_prob_all_roc <- do.call(bind_rows, exp_data_hap_prob_all) %>%
  mutate(HaplotypeProbability = signif(HaplotypeProbability, 2)) %>%
  group_by(HaplotypeProbability, Reads, Method, Pantranscriptome) %>%
  summarise(TP = sum(TP),
            TN = sum(TN),
            FP = sum(FP),
            FN = sum(FN),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  arrange(desc(HaplotypeProbability), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         TN_cs = cumsum(TN),
         FP_cs = cumsum(FP),
         FN_cs = cumsum(FN),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) %>%
  mutate(TPR_tpm = TP_cs / max(TP_cs + FN_cs), FPR_tpm = FP_cs / max(FP_cs + TN_cs)) %>%
  mutate(PPV_tpm = TP_cs / (TP_cs + FP_cs)) %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) 
  
exp_data_hap_prob_all_roc_points <- exp_data_hap_prob_all_roc %>% 
  mutate(HaplotypeProbability_floor = floor(HaplotypeProbability * 20) / 20) %>%
  group_by(HaplotypeProbability_floor, Reads, Method, Pantranscriptome) %>%
  mutate(is_min = (HaplotypeProbability == min(HaplotypeProbability))) %>%
  filter(is_min) %>%
  mutate(HaplotypeProbability = HaplotypeProbability_floor) %>%
  filter(HaplotypeProbability_floor == 0 | HaplotypeProbability_floor == 0.7 | HaplotypeProbability_floor == 0.8 | HaplotypeProbability_floor == 1)


exp_data_hap_exp_all_roc <- do.call(bind_rows, exp_data_hap_exp_all) %>%
  mutate(tpm_est_sig = signif(tpm_est, 1)) %>%
  mutate(tpm_est_sig = ifelse(tpm_est_sig >= 10, 10, tpm_est_sig)) %>%
  group_by(tpm_est_sig, Reads, Method, Pantranscriptome) %>%
  summarise(TP = sum(TP),
            TN = sum(TN),
            FP = sum(FP),
            FN = sum(FN),
            TP_tpm = sum(TP_tpm),
            FP_tpm = sum(FP_tpm)) %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  arrange(desc(tpm_est_sig), .by_group = T) %>%
  mutate(TP_cs = cumsum(TP),
         TN_cs = cumsum(TN),
         FP_cs = cumsum(FP),
         FN_cs = cumsum(FN),
         TP_tpm_cs = cumsum(TP_tpm),
         FP_tpm_cs = cumsum(FP_tpm)) %>%
  mutate(TPR_tpm = TP_cs / max(TP_cs + FN_cs), FPR_tpm = FP_cs / max(FP_cs + TN_cs)) %>%
  mutate(PPV_tpm = TP_cs / (TP_cs + FP_cs)) %>%
  mutate(frac_correct_tpm = TP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs), frac_incorrect_tpm = FP_tpm_cs / max(TP_tpm_cs + FP_tpm_cs)) %>%
  filter(tpm_est_sig > 0)

exp_data_hap_exp_all_roc_points <- exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  mutate(is_min = (tpm_est_sig == min(tpm_est_sig))) %>%
  mutate(is_max = (tpm_est_sig == max(tpm_est_sig))) %>%
  filter(is_min | tpm_est_sig == 0.1 | tpm_est_sig == 1 | is_max) 

exp_data_stats_all_bars <- do.call(bind_rows, exp_data_stats_all) %>%
  filter(Type != "Transcript") 


######


# wes_cols <- c(wes_palette("GrandBudapest1")[1], wes_palette("GrandBudapest2")[4], wes_palette("GrandBudapest1")[2], wes_palette("Chevalier1")[c(1,2)])
# 
# 
# exp_data_hap_prob_all_roc$FacetRow <- ""
# exp_data_hap_prob_all_roc_points$FacetRow <- ""
# 
# exp_data_hap_prob_all_roc_points_sim <- exp_data_hap_prob_all_roc_points %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/sim_hap_prob_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_prob_all_roc %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = HaplotypeProbability)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_prob_all_roc_points_sim, size = 2) +
#   geom_label_repel(data = exp_data_hap_prob_all_roc_points_sim, size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
#   scale_color_manual(values = wes_cols[c(4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("Expressed transcript precision") +
#   ylab("Expressed transcript sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# exp_data_hap_prob_all_roc_points_real <- exp_data_hap_prob_all_roc_points %>%
#   filter(Reads == "Real reads") %>%
#   filter(Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/real_hap_prob_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_prob_all_roc %>%
#   filter(Reads == "Real reads") %>%
#   filter(Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = HaplotypeProbability)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_prob_all_roc_points_real, size = 2) +
#   geom_label_repel(data = exp_data_hap_prob_all_roc_points_real, size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
#   scale_color_manual(values = wes_cols[c(4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   scale_x_continuous(trans = 'log10') +
#   scale_y_continuous(trans = 'log10') +
#   annotation_logticks() +
#   xlab("Number of expressed non-NA12878 transcripts") +
#   ylab("Number of expressed NA12878 transcripts") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# 
# exp_data_hap_exp_all_roc$FacetRow <- ""
# exp_data_hap_exp_all_roc_points$FacetRow <- ""
# 
# exp_data_hap_exp_all_roc_points_sim <- exp_data_hap_exp_all_roc_points %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/sim_expression_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_exp_all_roc %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_exp_all_roc_points_sim, size = 2) +
#   geom_label_repel(data = exp_data_hap_exp_all_roc_points_sim[exp_data_hap_exp_all_roc_points_sim$Pantranscriptome == "Whole (excl. CEU)" & exp_data_hap_exp_all_roc_points_sim$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
#   scale_color_manual(values = wes_cols[c(1,2,3,4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("Transcript expression precision") +
#   ylab("Transcript expression sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# exp_data_hap_exp_all_roc_points_sim_sample <- exp_data_hap_exp_all_roc_points %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome == "Sample-specific (NA12878)")
# 
# pdf(paste("plots/sim_expression_roc_pre_sample_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_exp_all_roc %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome == "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_exp_all_roc_points_sim_sample, size = 2) +
#   geom_label_repel(data = exp_data_hap_exp_all_roc_points_sim_sample, size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
#   scale_color_manual(values = wes_cols[c(1,2,3,4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("Transcript expression precision") +
#   ylab("Transcript expression sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# exp_data_hap_exp_all_roc_points_sim <- exp_data_hap_exp_all_roc_points %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "rpvg" | Method == "rpvg (single-path)") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/sim_expression_roc_pre_rpvg_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_exp_all_roc %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Method == "rpvg" | Method == "rpvg (single-path)") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_exp_all_roc_points_sim, size = 2) +
#   geom_label_repel(data = exp_data_hap_exp_all_roc_points_sim[exp_data_hap_exp_all_roc_points_sim$Pantranscriptome == "Whole (excl. CEU)" & exp_data_hap_exp_all_roc_points_sim$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
#   scale_color_manual(values = wes_cols[c(4,5)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("Transcript expression precision") +
#   ylab("Transcript expression sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# exp_data_hap_exp_all_roc_points_real <- exp_data_hap_exp_all_roc_points %>%
#   filter(Reads == "Real reads") %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/real_expression_roc_pre_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_exp_all_roc %>%
#   filter(Reads == "Real reads") %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig, segment.color = Method)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_exp_all_roc_points_real, size = 2) +
#   geom_label_repel(data = exp_data_hap_exp_all_roc_points_real[exp_data_hap_exp_all_roc_points_real$Pantranscriptome == "Whole (excl. CEU)" & exp_data_hap_exp_all_roc_points_real$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
#   scale_color_manual(values = wes_cols[c(1,2,3,4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   scale_x_continuous(trans = 'log10') +
#   scale_y_continuous(trans = 'log10') +
#   annotation_logticks() +
#   xlab("Number of expressed non-NA12878 transcripts") +
#   ylab("Number of expressed NA12878 transcripts") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# exp_data_hap_exp_all_roc_points_real <- exp_data_hap_exp_all_roc_points %>%
#   filter(Reads == "Real reads") %>%
#   filter(Method == "rpvg" | Method == "rpvg (single-path)") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/real_expression_roc_pre_rpvg_", dataset, ".pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_exp_all_roc %>%
#   filter(Reads == "Real reads") %>%
#   filter(Method == "rpvg" | Method == "rpvg (single-path)") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TP_cs, x = FP_cs, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_exp_all_roc_points_real, size = 2) +
#   geom_label_repel(data = exp_data_hap_exp_all_roc_points_real[exp_data_hap_exp_all_roc_points_real$Pantranscriptome == "Whole (excl. CEU)" & exp_data_hap_exp_all_roc_points_real$tpm_est_sig %in% c(0.1,1,10),], size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30") +
#   scale_color_manual(values = wes_cols[c(4,5)]) +
#   facet_grid(FacetRow ~ Reads) +
#   scale_x_continuous(trans = 'log10') +
#   scale_y_continuous(trans = 'log10') +
#   annotation_logticks() +
#   xlab("Number of expressed non-NA12878 transcripts") +
#   ylab("Number of expressed NA12878 transcripts") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# 
# exp_data_stats_all_bars_main <- exp_data_stats_all_bars %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   ungroup() %>%
#   add_row(Reads = "Simulated reads", Method = "RSEM", Pantranscriptome = "Whole\n(excl. CEU)", Type = "All", frac_hap_error_tpm = 1, Truncated = FALSE) %>%
#   add_row(Reads = "Simulated reads", Method = "RSEM", Pantranscriptome = "Whole", Type = "All", frac_hap_error_tpm = 1, Truncated = FALSE) %>%
#   add_row(Reads = "Real reads", Method = "RSEM", Pantranscriptome = "Whole\n(excl. CEU)", Type = "All", frac_hap_error_tpm = 1, Truncated = FALSE) %>%
#   add_row(Reads = "Real reads", Method = "RSEM", Pantranscriptome = "Whole", Type = "All", frac_hap_error_tpm = 1, Truncated = FALSE)
# 
# exp_data_stats_all_bars_main$Reads <- factor(exp_data_stats_all_bars_main$Reads, levels = c("Simulated reads", "Real reads"))
# exp_data_stats_all_bars_main$Method <- factor(exp_data_stats_all_bars_main$Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg", "rpvg (single-path)"))
# 
# exp_data_stats_all_bars_main$Pantranscriptome = recode_factor(exp_data_stats_all_bars_main$Pantranscriptome,
#                            "Sample-specific (NA12878)" = "Sample-specific\n(NA12878)",
#                            "Europe (excl. CEU)" = "Europe\n(excl. CEU)",
#                            "Whole (excl. CEU)" = "Whole\n(excl. CEU)",
#                            "Whole" = "Whole"
# )
# 
# exp_data_stats_all_bars_main$Pantranscriptome <- factor(exp_data_stats_all_bars_main$Pantranscriptome, levels = c("Sample-specific\n(NA12878)", "Europe\n(excl. CEU)", "Whole\n(excl. CEU)", "Whole"))
# 
# exp_data_stats_all_bars_main$FacetRow <- ""
# 
# pdf(paste("plots/sim_real_hap_tpm_error_", dataset, ".pdf", sep = ""), height = 4, width = 7, pointsize = 12)
# exp_data_stats_all_bars_main %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg") %>%
#   filter(Pantranscriptome != "Sample-specific\n(NA12878)") %>%
#   filter(Type == "All") %>%
#   filter(!Truncated) %>%
#   ggplot(aes(x = Pantranscriptome, y = 1 - frac_hap_error_tpm, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
#   scale_fill_manual(values = wes_cols[c(1,2,3,4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("") +
#   ylab("Fraction TPM on NA12878 haplotypes") +
#   scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 12))
# dev.off()
# 
# pdf(paste("plots/sim_real_hap_tpm_error_rpvg_", dataset, ".pdf", sep = ""), height = 4, width = 7.5, pointsize = 12)
# exp_data_stats_all_bars_main %>%
#   filter(Method == "rpvg" | Method == "rpvg (single-path)") %>%
#   filter(Pantranscriptome != "Sample-specific\n(NA12878)") %>%
#   filter(Type == "All") %>%
#   filter(!Truncated) %>%
#   ggplot(aes(x = Pantranscriptome, y = 1 - frac_hap_error_tpm, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
#   scale_fill_manual(values = wes_cols[c(4,5)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("") +
#   ylab("Fraction TPM on NA12878 haplotypes") +
#   scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 12))
# dev.off()
# 
# exp_data_stats_all_bars_main_sim <- exp_data_stats_all_bars_main %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Type != "Transcript") %>% 
#   filter(!Truncated)
# 
# exp_data_stats_all_bars_main_sim <- exp_data_stats_all_bars_main_sim %>%
#    ungroup() %>%
#    add_row(Reads = "Simulated reads", Method = "RSEM", Pantranscriptome = "Whole\n(excl. CEU)", Type = "All", Spearman_tpm = 0, ARD_mean_tpm = 0, Truncated = FALSE) %>%
#    add_row(Reads = "Simulated reads", Method = "RSEM", Pantranscriptome = "Whole\n(excl. CEU)", Type = "Haplotype", Spearman_tpm = 0, ARD_mean_tpm = 0, Truncated = FALSE) %>%
#    add_row(Reads = "Simulated reads", Method = "RSEM", Pantranscriptome = "Whole", Type = "All", Spearman_tpm = 0, ARD_mean_tpm = 0, Truncated = FALSE) %>%
#    add_row(Reads = "Simulated reads", Method = "RSEM", Pantranscriptome = "Whole", Type = "Haplotype", Spearman_tpm = 0, ARD_mean_tpm = 0, Truncated = FALSE)
# 
# exp_data_stats_all_bars_main_sim$Type <- as.factor(exp_data_stats_all_bars_main_sim$Type)
# exp_data_stats_all_bars_main_sim$Method <- factor(exp_data_stats_all_bars_main_sim$Method, levels = c("Kallisto", "Salmon", "RSEM", "rpvg", "rpvg (single-path)"))
# exp_data_stats_all_bars_main_sim$Pantranscriptome <- factor(exp_data_stats_all_bars_main_sim$Pantranscriptome, levels = c("Sample-specific\n(NA12878)", "Europe\n(excl. CEU)", "Whole\n(excl. CEU)", "Whole"))
# 
# exp_data_stats_all_bars_main_sim$FacetRow <- ""
# 
# exp_data_stats_all_bars_main_sim_nosingle <- exp_data_stats_all_bars_main_sim %>%
#   filter(Method == "Kallisto" | Method == "Salmon" | Method == "RSEM" | Method == "rpvg")
#   
# pdf(paste("plots/sim_mard_all_hap_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
# ggplot() +
#   geom_bar(data = subset(exp_data_stats_all_bars_main_sim_nosingle, Type == "All"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
#   geom_bar(data = subset(exp_data_stats_all_bars_main_sim_nosingle, Type == "Haplotype"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.6) +
#   scale_alpha_manual(name = "Transcripts", values = c(1, 0.6), labels = c("All", "NA12878"), drop = F) +
#   scale_fill_manual(values = wes_cols[c(1,2,3,4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("") +
#   ylab("Mean relative expression difference") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 13))
# dev.off()
# 
# pdf(paste("plots/sim_spearman_all_hap_tpm_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
# ggplot() +
#   geom_bar(data = subset(exp_data_stats_all_bars_main_sim_nosingle, Type == "All"), aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
#   geom_bar(data = subset(exp_data_stats_all_bars_main_sim_nosingle, Type == "Haplotype"), aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.6) +
#   scale_alpha_manual(name = "Transcripts", values = c(1, 0.6), labels = c("All", "NA12878"), drop = F) +
#   scale_fill_manual(values = wes_cols[c(1,2,3,4)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("") +
#   ylab("Spearman expression correlation") +
#   scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 13))
# dev.off()
# 
# exp_data_stats_all_bars_main_sim_rpvg <- exp_data_stats_all_bars_main_sim %>%
#   filter(Method == "rpvg" | Method == "rpvg (single-path)")
# 
# pdf(paste("plots/sim_mard_all_hap_tpm_rpvg_", dataset, ".pdf", sep = ""), height = 4, width = 6, pointsize = 12)
# ggplot() +
#   geom_bar(data = subset(exp_data_stats_all_bars_main_sim_rpvg, Type == "All"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge()) +
#   geom_bar(data = subset(exp_data_stats_all_bars_main_sim_rpvg, Type == "Haplotype"), aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method, alpha = Type), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.6) +
#   scale_alpha_manual(name = "Transcripts", values = c(1, 0.6), labels = c("All", "NA12878"), drop = F) +
#   scale_fill_manual(values = wes_cols[c(4,5)]) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("") +
#   ylab("Mean relative expression difference") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 13))
# dev.off()
# 
# }

##### Debug 


exp_data_stats_all_bars %>%
  filter(Pantranscriptome != "NA12878") %>%
  filter(Type == "All") %>%
  select(Reads, Method, Pantranscriptome, n, frac_hap_error_count, frac_hap_error_tpm) %>%
  print(n = 100)

exp_data_hap_exp_all_roc %>%
  group_by(Reads, Method, Pantranscriptome) %>%
  filter(tpm_est_sig == 0.01) %>%
  summarise(max_FP_cs = max(FP_cs)) %>%
  print(n = 100)


pdf(paste("plots/hap_prob_rocs_pt_debug_", dataset, ".pdf", sep = ""), width = 9, height = 9)

for (pt in unique(exp_data_hap_prob_all_roc$Pantranscriptome)) {
  
  exp_data_hap_prob_all_roc_points_pt <- exp_data_hap_prob_all_roc_points %>%
    filter(Pantranscriptome == pt)

  delta <- max(max(exp_data_hap_prob_all_roc$PPV_tpm) - min(exp_data_hap_prob_all_roc$PPV_tpm), max(exp_data_hap_prob_all_roc$TPR_tpm) - min(exp_data_hap_prob_all_roc$TPR_tpm))
  
  p <- exp_data_hap_prob_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$PPV_tpm), min(exp_data_hap_prob_all_roc$PPV_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$TPR_tpm), min(exp_data_hap_prob_all_roc$TPR_tpm) + delta)) +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression sensitivity (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  delta <- max(max(exp_data_hap_prob_all_roc$FP_cs) - min(exp_data_hap_prob_all_roc$FP_cs), max(exp_data_hap_prob_all_roc$TP_cs) - min(exp_data_hap_prob_all_roc$TP_cs))
  
  p <- exp_data_hap_prob_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Method, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads) +
    scale_x_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$FP_cs), min(exp_data_hap_prob_all_roc$FP_cs) + delta)) +
    scale_y_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$TP_cs), min(exp_data_hap_prob_all_roc$TP_cs) + delta)) +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  delta <- max(max(exp_data_hap_prob_all_roc$frac_incorrect_tpm) - min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), max(exp_data_hap_prob_all_roc$frac_correct_tpm) - min(exp_data_hap_prob_all_roc$frac_correct_tpm))
  
  p <- exp_data_hap_prob_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Method, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), min(exp_data_hap_prob_all_roc$frac_incorrect_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$frac_correct_tpm), min(exp_data_hap_prob_all_roc$frac_correct_tpm) + delta)) +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/hap_prob_rocs_me_debug_", dataset, ".pdf", sep = ""), width = 9, height = 9)

for (me in unique(exp_data_hap_prob_all_roc$Method)) {
  
  exp_data_hap_prob_all_roc_points_me <- exp_data_hap_prob_all_roc_points %>%
    filter(Method == me)
  
  delta <- max(max(exp_data_hap_prob_all_roc$PPV_tpm) - min(exp_data_hap_prob_all_roc$PPV_tpm), max(exp_data_hap_prob_all_roc$TPR_tpm) - min(exp_data_hap_prob_all_roc$TPR_tpm))
  
  p <- exp_data_hap_prob_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Pantranscriptome, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$PPV_tpm), min(exp_data_hap_prob_all_roc$PPV_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$TPR_tpm), min(exp_data_hap_prob_all_roc$TPR_tpm) + delta)) +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression sensitivity (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  delta <- max(max(exp_data_hap_prob_all_roc$FP_cs) - min(exp_data_hap_prob_all_roc$FP_cs), max(exp_data_hap_prob_all_roc$TP_cs) - min(exp_data_hap_prob_all_roc$TP_cs))
  
  p <- exp_data_hap_prob_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Pantranscriptome, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads) +
    scale_x_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$FP_cs), min(exp_data_hap_prob_all_roc$FP_cs) + delta)) +
    scale_y_continuous(trans = 'log10', limits = c(min(exp_data_hap_prob_all_roc$TP_cs), min(exp_data_hap_prob_all_roc$TP_cs) + delta)) +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  delta <- max(max(exp_data_hap_prob_all_roc$frac_incorrect_tpm) - min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), max(exp_data_hap_prob_all_roc$frac_correct_tpm) - min(exp_data_hap_prob_all_roc$frac_correct_tpm))
  
  p <- exp_data_hap_prob_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Pantranscriptome, label = HaplotypeProbability)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_prob_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_prob_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads) +
    xlim(c(min(exp_data_hap_prob_all_roc$frac_incorrect_tpm), min(exp_data_hap_prob_all_roc$frac_incorrect_tpm) + delta)) +
    ylim(c(min(exp_data_hap_prob_all_roc$frac_correct_tpm), min(exp_data_hap_prob_all_roc$frac_correct_tpm) + delta)) +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/expression_rocs_pt_debug_", dataset, ".pdf", sep = ""))

for (pt in unique(exp_data_hap_exp_all_roc$Pantranscriptome)) {

  exp_data_hap_exp_all_roc_points_pt <- exp_data_hap_exp_all_roc_points %>%
    filter(Pantranscriptome == pt)
  
  p <- exp_data_hap_exp_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads, scales="free") +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression sensitivity (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  p <- exp_data_hap_exp_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Method, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads, scales="free") +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  p <- exp_data_hap_exp_all_roc %>%
    filter(Pantranscriptome == pt) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Method, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_pt, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_pt, size = 2, fontface = 2) +
    facet_wrap(Pantranscriptome~Reads, scales="free") +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
#    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/expression_rocs_me_debug_", dataset, ".pdf", sep = ""))

for (me in unique(exp_data_hap_exp_all_roc$Method)) {
  
  exp_data_hap_exp_all_roc_points_me <- exp_data_hap_exp_all_roc_points %>%
    filter(Method == me)
  
  p <- exp_data_hap_exp_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Pantranscriptome, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads, scales="free") +
    xlab("Transcript expression precision (PPV)") +
    ylab("Transcript expression sensitivity (TPR)") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  p <- exp_data_hap_exp_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = TP_cs, x = FP_cs, color = Pantranscriptome, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads, scales="free") +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab("Number incorrect expressed transcripts") +
    ylab("Number correct expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
  
  p <- exp_data_hap_exp_all_roc %>%
    filter(Method == me) %>%
    ggplot(aes(y = frac_correct_tpm, x = frac_incorrect_tpm, color = Pantranscriptome, label = tpm_est_sig)) +
    geom_line(size = 0.5) +
    geom_point(data = exp_data_hap_exp_all_roc_points_me, size = 1.5) +
    geom_text_repel(data = exp_data_hap_exp_all_roc_points_me, size = 2, fontface = 2) +
    facet_wrap(Method~Reads, scales="free") +
    xlab("Fraction TPM on not expressed transcripts") +
    ylab("Fraction TPM on expressed transcripts") +
    guides(linetype = FALSE) +
    guides(shape = FALSE) +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    #    theme(legend.position="bottom") +
    theme(text = element_text(size=12))
  print(p)
}

dev.off()


pdf(paste("plots/general_stats_debug_", dataset, ".pdf", sep = ""))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = frac_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Fraction estimated TPM for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = frac_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Fraction estimated CPM for transcripts\non incorrect haplotypes") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = num_hap_error_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = num_hap_error_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  scale_y_continuous(trans = 'log10') +
  annotation_logticks() +
  xlab("") +
  ylab("Number of transcripts on incorrect\nhaplotypes (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Type == "All") %>%
  ggplot(aes(x = Pantranscriptome, y = ExpCorrect, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Fraction correct expressed") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads") %>%
  ggplot(aes(x = Pantranscriptome, y = Pearson_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads") %>%
  ggplot(aes(x = Pantranscriptome, y = Pearson_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Pearson correlation (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads") %>%
  ggplot(aes(x = Pantranscriptome, y = Spearman_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads") %>%
  ggplot(aes(x = Pantranscriptome, y = Spearman_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  #  scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
  xlab("") +
  ylab("Spearman correlation (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads") %>%
  ggplot(aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Mean absolute difference (TPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

exp_data_stats_all_bars %>%
  filter(Reads == "Simulated reads") %>%
  ggplot(aes(x = Pantranscriptome, y = ARD_mean_count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
  facet_grid(Reads ~ Type, scales="free") +
  xlab("") +
  ylab("Mean absolute difference (CPM)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=8))

dev.off()




##### Salmon evaluation
# 
# exp_data_hap_exp_all_roc_points_sim_salmon <- exp_data_hap_exp_all_roc_points %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)")
# 
# pdf(paste("plots/sim_expression_roc_pre_", dataset, "_salmon.pdf", sep = ""), height = 5, width = 7, pointsize = 12)
# exp_data_hap_exp_all_roc %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Pantranscriptome != "Sample-specific (NA12878)") %>%
#   ggplot(aes(y = TPR_tpm, x = PPV_tpm, color = Method, linetype = Pantranscriptome, shape = Pantranscriptome, label = tpm_est_sig)) +
#   geom_line(size = 1) +
#   geom_point(data = exp_data_hap_exp_all_roc_points_sim_salmon, size = 2) +
#   geom_label_repel(data = exp_data_hap_exp_all_roc_points_sim_salmon, size = 2.5, fontface = 2, box.padding = 0.5, min.segment.length = 0, show.legend = FALSE, segment.color = "grey30", xlim = c(0, 0.99)) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("Transcript expression precision") +
#   ylab("Transcript expression sensitivity") +
#   theme_bw() +
#   theme(aspect.ratio = 1) +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(legend.key.width = unit(1, "cm")) +
#   theme(text = element_text(size = 15))
# dev.off()
# 
# exp_data_stats_all_bars_salmon <- exp_data_stats_all_bars %>%
#   add_column(FacetRow = "") %>%
#   filter(Reads == "Simulated reads") %>%
#   filter(Type != "Transcript") 
# 
# pdf(paste("plots/sim_real_hap_tpm_error_", dataset, "_salmon.pdf", sep = ""), height = 4, width = 7, pointsize = 12)
# exp_data_stats_all_bars_salmon %>%
#   filter(Type == "All") %>%
#   ggplot(aes(x = Pantranscriptome, y = 1 - frac_hap_error_tpm, fill = Method)) +
#   geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
#   facet_grid(FacetRow ~ Reads) +
#   xlab("") +
#   ylab("Fraction TPM on NA12878 haplotypes") +
#   scale_y_continuous(limits = c(0, 1), oob = rescale_none) +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 12))
# dev.off()
# 
# exp_data_stats_all_bars_salmon$Type = recode_factor(exp_data_stats_all_bars_salmon$Type,
#                             "Haplotype" = "NA12878")
# exp_data_stats_all_bars_salmon$Type <- factor(exp_data_stats_all_bars_salmon$Type, levels = c("All", "NA12878"))
# 
# pdf(paste("plots/sim_mard_all_hap_tpm_", dataset, "_salmon.pdf", sep = ""), height = 4, width = 6, pointsize = 12)
# exp_data_stats_all_bars_salmon %>%
#   ggplot() +
#   geom_bar(aes(x = Pantranscriptome, y = ARD_mean_tpm, fill = Method), stat = "identity", width = 0.5, position = position_dodge()) +
#   facet_grid(FacetRow ~ Type) +
#   xlab("") +
#   ylab("Mean relative expression difference") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(panel.spacing = unit(0.5, "cm")) +
#   theme(text = element_text(size = 13))
# dev.off()
