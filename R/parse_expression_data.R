
# parse_expression_data.R

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

parse_salmon <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypePosterior = 1) %>%
    rename(name = Name, tpm_est = TPM, count_est = NumReads) %>%
    select(-Length, -EffectiveLength)
  
  return(data)
}

parse_kallisto <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypePosterior = 1) %>%
    rename(name = target_id, tpm_est = tpm, count_est = est_counts) %>%
    select(-length, -eff_length)
  
  return(data)
}

parse_rsem <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypePosterior = 1) %>%
    rename(name = transcript_id, tpm_est = TPM, count_est = expected_count) %>%
    select(-gene_id, -length, -effective_length, -FPKM, -IsoPct)
  
  return(data)
}

parse_rpvg <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    rename(name = Name, tpm_est = TPM, count_est = ReadCount) %>%
    select(-Length, -EffectiveLength, -ClusterRelativeExpression)
  
  if (grepl("_multi_", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "multi", sep = "_"))
  } 
  
  return(data)
}

getStats <- function(data) {

  data_stats <- data %>%
    mutate(ard_count = abs(count_est - count_sim) / (count_est + count_sim)) %>%
    mutate(ard_tpm = abs(tpm_est - tpm_sim) / (tpm_est + tpm_sim)) %>%
    replace_na(list(ard_count = 0)) %>%
    replace_na(list(ard_tpm = 0)) %>%
    group_by(Reads, Method, Graph) %>%
    summarise(
      n = n(), 
      Expressed = sum(tpm_est > 0), 
      ExpCorrect = mean((tpm_est > 0) == (tpm_sim > 0)), 
      ExpSimCorrect = sum((tpm_est > 0) & (tpm_sim > 0)) / sum(tpm_sim > 0), 
      Pearson_count = cor(count_est, count_sim, method = "pearson"),
      Pearson_tpm = cor(tpm_est, tpm_sim, method = "pearson"),
      Spearman_count = cor(count_est, count_sim, method = "spearman"),
      Spearman_tpm = cor(tpm_est, tpm_sim, method = "spearman"),
      ARD_mean_count = mean(ard_count),
      ARD_mean_tpm = mean(ard_tpm),
      frac_hap_error_count = sum((!is_hap) * count_est) / sum(count_est),
      frac_hap_error_tpm = sum((!is_hap) * tpm_est) / sum(tpm_est),
      num_hap_error_count = sum((!is_hap) & count_est > 0),
      num_hap_error_tpm = sum((!is_hap) & tpm_est > 0))
  
  return(data_stats)
}

# sim_exp_h1 <- read_table2("sim/SRR1153470/vg/sim_1kg_NA12878_gencode100_SRR1153470_vg_h1_paths.txt.gz") %>%
#   group_by(path) %>%
#   summarise(count = n()) %>%
#   mutate(count = count / 2)
# 
# sim_exp_h2 <- read_table2("sim/SRR1153470/vg/sim_1kg_NA12878_gencode100_SRR1153470_vg_h2_paths.txt.gz") %>%
#   group_by(path) %>%
#   summarise(count = n()) %>%
#   mutate(count = count / 2)
# 
# sim_exp <- bind_rows(sim_exp_h1, sim_exp_h2)
# 
# save(sim_exp, file = "sim/SRR1153470/vg/sim_1kg_NA12878_gencode100_SRR1153470_vg.RData")

read_type <- "sim_vg"

identical_seqs <- read_table2("graphs/1kg_NA12878_exons_gencode100_allpaths/1kg_NA12878_gencode100_genes_hst_overlap.txt")
#identical_seqs <- read_table2("graphs/1kg_NA12878_exons_gencode100_allpaths/1kg_nonCEU_af001_gencode100_genes_hst_overlap.txt")
#identical_seqs <- read_table2("graphs/1kg_NA12878_exons_gencode100_allpaths/1kg_all_af001_gencode100_genes_hst_overlap.txt")

rsem <- read_table2("rsem/SRR1153470/1kg_NA12878_gencode100_SRR1153470_rsem.isoforms.results")

load("sim/SRR1153470/vg/sim_1kg_NA12878_gencode100_SRR1153470_vg.RData")

sim_exp <- sim_exp %>%
  right_join(rsem, by = c("path" = "transcript_id")) %>%
  left_join(identical_seqs, by = c("path" = "Name1")) %>%
  rename(name = Name2) %>%
  mutate(name = ifelse(is.na(name), paste(path, "_na", sep = ""), name)) %>%
  mutate(TPM = count / effective_length) %>%
  replace_na(list(TPM = 0, count = 0)) %>%
  mutate(TPM = 10^6 * TPM / sum(TPM)) %>%
  group_by(name) %>%
  summarise(length = length, effective_length = max(effective_length), tpm_sim = sum(TPM), count_sim = sum(count)) 

sim_exp %>% filter(tpm_sim > 0) %>% arrange(tpm_sim)

sim_exp %>% filter(tpm_sim > 0 & count_sim == 0) %>% print()
sim_exp %>% filter(tpm_sim == 0 & count_sim > 0) %>% print()

if (read_type == "real") {
  
  sim_exp <- sim_exp %>%
    mutate(tpm_sim = 1) %>%
    mutate(count_sim = 1)
}

for (f in list.files(path = "methods", pattern = "rpvg9_exact_mpmap3_multi.*.gz", full.names = T, recursive = T)) { 

  if (!grepl(read_type, f)) {
    
    next 
  }
    
  if (!grepl("1kg_NA12878_gencode100", f)) {
    
    next 
  }
  
  if (grepl("quant.sf", f)) {
    
    exp_data <- parse_salmon(f)
  
  } else if (grepl("abundance.tsv", f)) {

    exp_data <- parse_kallisto(f)

  } else if (grepl("isoforms.results", f)) {

    exp_data <- parse_rsem(f)
    
  } else {
    
    exp_data <- parse_rpvg(f)
  }

  print(f)
  
  exp_data <- exp_data %>% 
    full_join(sim_exp, by = "name") %>%
    mutate(is_hap = ifelse(is.na(tpm_sim), F, T)) %>%
    replace_na(list(HaplotypePosterior = 0, tpm_est = 0, count_est = 0, Reads = exp_data$Reads[1], Method = exp_data$Method[1], Graph = exp_data$Graph[1], tpm_sim = 0, count_sim = 0)) %>%
    separate(name, c("transcript", "hap_id"), "_")
  
  exp_data %>% filter(tpm_est > 0 & count_est == 0) %>% print()
  exp_data %>% filter(tpm_est == 0 & count_est > 0) %>% print()
  
  # exp_data %>% filter(count_sim > 0 & count_est == 0) %>% arrange(desc(count_sim))
  # exp_data %>% arrange(desc(abs(count_sim - count_est)))
  # 
  # exp_data %>% filter(count_est > 0) %>% ungroup() %>% summarise(min_count_est = min(count_est), min_tpm_est = min(tpm_est)) %>% print()
  # 
  
  exp_data %>% filter(!is_hap) %>% filter(HaplotypePosterior >= 0.5) %>% filter(tpm_est > 0) %>% arrange(desc(tpm_est)) %>% mutate(rel_tpm_est = tpm_est / sum(tpm_est))

  exp_data %>% filter(transcript == "ENST00000646664.1") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000227378.7") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000514057.1") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000394667.7") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000253788.11") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000511473.5") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000287038.7") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000370321.8") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000412585.6") %>% filter((!is_hap & count_est > 0) | is_hap)
  exp_data %>% filter(transcript == "ENST00000648437.1") %>% filter((!is_hap & count_est > 0) | is_hap)
  
  
  
  print(nrow(exp_data))
  
  exp_data_hap_pos <- exp_data %>%
    group_by(HaplotypePosterior, Reads, Method, Graph) %>%
    summarise(TP = sum((tpm_sim > 0) & (tpm_est > 0)),
              TN = sum((tpm_sim == 0) & (tpm_est == 0)),
              FP = sum((tpm_sim == 0) & (tpm_est > 0)),
              FN = sum((tpm_sim > 0) & (tpm_est == 0)),
              TP_count = sum((tpm_sim > 0) * count_est),
              FP_count = sum((tpm_sim == 0) * count_est),
              TP_tpm = sum((tpm_sim > 0) * tpm_est),
              FP_tpm = sum((tpm_sim == 0) * tpm_est))
  
  print(nrow(exp_data_hap_pos))
  
  exp_data <- exp_data %>% 
    mutate(count_est = ifelse(HaplotypePosterior >= 0.5, count_est, 0)) %>%
    mutate(tpm_est = ifelse(HaplotypePosterior >= 0.5, tpm_est, 0))
  
  exp_data_hap_exp <- exp_data %>%
    group_by(tpm_est, Reads, Method, Graph) %>%
    summarise(TP = sum((tpm_sim > 0) & (tpm_est > 0)),
              TN = sum((tpm_sim == 0) & (tpm_est == 0)),
              FP = sum((tpm_sim == 0) & (tpm_est > 0)),
              FN = sum((tpm_sim > 0) & (tpm_est == 0)),
              TP_count = sum((tpm_sim > 0) * count_est),
              FP_count = sum((tpm_sim == 0) * count_est),
              TP_tpm = sum((tpm_sim > 0) * tpm_est),
              FP_tpm = sum((tpm_sim == 0) * tpm_est))
  
  print(nrow(exp_data_hap_exp))
  
  exp_data_stats_all <- exp_data %>%
    getStats() %>%
    add_column(Type = "All")

  exp_data_stats_hap <- exp_data %>%
    filter(is_hap) %>%
    getStats() %>%
    add_column(Type = "Haplotype")
  
  exp_data_stats_txp <- exp_data %>%
    group_by(transcript, Reads, Method, Graph) %>%
    summarise(count_est = sum(count_est), count_sim = sum(count_sim), tpm_est = sum(tpm_est), tpm_sim = sum(tpm_sim), is_hap = max(is_hap)) %>%
    ungroup() %>%
    getStats() %>%
    add_column(Type = "Transcript")
  
  exp_data_stats <- rbind(exp_data_stats_all, exp_data_stats_hap, exp_data_stats_txp) %>%
    add_column(Truncated = FALSE)

  exp_data_trunc <- exp_data %>%
    mutate(count_est = ifelse(count_est >= 0.1, count_est, 0)) %>%  
    mutate(tpm_est = ifelse(tpm_est >= 10^-4, tpm_est, 0))
    
  exp_data_stats_all_trunc <- exp_data_trunc %>%
    getStats() %>%
    add_column(Type = "All")
  
  exp_data_stats_hap_trunc <- exp_data_trunc %>%
    filter(is_hap) %>%
    getStats() %>%
    add_column(Type = "Haplotype")
  
  exp_data_stats_txp_trunc <- exp_data_trunc %>%
    group_by(transcript, Reads, Method, Graph) %>%
    summarise(count_est = sum(count_est), count_sim = sum(count_sim), tpm_est = sum(tpm_est), tpm_sim = sum(tpm_sim), is_hap = max(is_hap)) %>%
    ungroup() %>%
    getStats() %>%
    add_column(Type = "Transcript")
  
  exp_data_stats_trunc <- rbind(exp_data_stats_all_trunc, exp_data_stats_hap_trunc, exp_data_stats_txp_trunc) %>%
    add_column(Truncated = TRUE)
  
  exp_data_stats <- rbind(exp_data_stats, exp_data_stats_trunc)
  print(exp_data_stats)

  save(exp_data_hap_pos, exp_data_hap_exp, exp_data_stats, file = paste("rdata/", exp_data_stats$Method[1], exp_data_stats$Reads[1], exp_data_stats$Graph[1], ".RData", sep = "", collapse = ""))
}
