
# parse_expression_data.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")
library("truncnorm")

source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/quant_debug/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########


hap_prob_thres <- 0.8

update_tpm <- function(exp_data) {
  
  exp_data <- exp_data %>%
    mutate(TPM = count / (length - etruncnorm(1, length, sim_mean, sim_sd))) %>%
    replace_na(list(TPM = 0)) %>%
    mutate(TPM = 10^6 * TPM / sum(TPM))
  
  return(exp_data)
}

parse_salmon <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = Name, tpm_est = TPM, count_est = NumReads, length = Length) %>%
    select(-EffectiveLength)
  
  if ("FracNonZero" %in% names(data)) {
    
    data <- data %>%
      select(-FracNonZero)     
  }
  
  return(data)
}

parse_kallisto <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = target_id, tpm_est = tpm, count_est = est_counts) %>%
    select(-eff_length)
  
  return(data)
}

parse_rsem <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    add_column(HaplotypeProbability = 1) %>%
    rename(name = transcript_id, tpm_est = TPM, count_est = expected_count) %>%
    select(-gene_id, -effective_length, -FPKM, -IsoPct)
  
  if ("posterior_mean_count" %in% names(data)) {
    
    data <- data %>%
      mutate(count_est = posterior_mean_count) %>%
      mutate(tpm_est = pme_TPM) %>%
      select(-posterior_mean_count, -posterior_standard_deviation_of_count, -pme_TPM, -pme_FPKM, -IsoPct_from_pme_TPM)     
  }
  
  return(data)
}

parse_rpvg <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(gzfile(filename), col_names = T)
  data <- data %>%
    add_column(Reads = paste(dir_split[5], dir_split[6], sep = "_")) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8]) %>%
    rename(name = Name, tpm_est = TPM, count_est = ReadCount, length = Length) %>%
    select(-ClusterID, -EffectiveLength)
  
  if (grepl("map_fast", basename(filename))) {
    
    data <- data %>%
      mutate(Method = paste(Method, "map", sep = "_"))
  }
  
  return(data)
}


getStats <- function(data) {

  data_stats <- data %>%
    mutate(ard_count = abs(count_est - count_truth) / (count_est + count_truth)) %>%
    mutate(ard_tpm = abs(tpm_est - tpm_truth) / (tpm_est + tpm_truth)) %>%
    replace_na(list(ard_count = 0)) %>%
    replace_na(list(ard_tpm = 0)) %>%
    group_by(Reads, Method, Graph) %>%
    summarise(
      n = n(), 
      Expressed = sum(tpm_est > 0), 
      ExpCorrect = mean((tpm_est > 0) == (tpm_truth > 0)), 
      ExpSimCorrect = sum((tpm_est > 0) & (tpm_truth > 0)) / sum(tpm_truth > 0), 
      Pearson_count = cor(count_est, count_truth, method = "pearson"),
      Pearson_tpm = cor(tpm_est, tpm_truth, method = "pearson"),
      Spearman_count = cor(count_est, count_truth, method = "spearman"),
      Spearman_tpm = cor(tpm_est, tpm_truth, method = "spearman"),
      ARD_mean_count = mean(ard_count),
      ARD_mean_tpm = mean(ard_tpm),
      frac_hap_error_count = sum((!is_hap) * count_est) / sum(count_est),
      frac_hap_error_tpm = sum((!is_hap) * tpm_est) / sum(tpm_est),
      num_hap_error_count = sum((!is_hap) & count_est > 0),
      num_hap_error_tpm = sum((!is_hap) & tpm_est > 0))
  
  return(data_stats)
}

parse_data <- function(dataset, sim_mean, sim_sd, read_type, ref_name) {
 
  # truth_exp_h1 <- read_table2(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1_h1_paths.txt.gz", sep = "")) %>%
  #   group_by(path) %>%
  #   summarise(count = n()) %>%
  #   mutate(count = count / 2)
  # 
  # truth_exp_h2 <- read_table2(paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1_h2_paths.txt.gz", sep = "")) %>%
  #   group_by(path) %>%
  #   summarise(count = n()) %>%
  #   mutate(count = count / 2)
  # 
  # truth_exp <- bind_rows(truth_exp_h1, truth_exp_h2)
  # save(truth_exp, file = paste("sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1.RData", sep = ""))
  
  identical_seqs <- read_table2(paste("../quant_r1/graphs/1kg_NA12878_exons_gencode100_allpaths/", ref_name, "_hst_overlap.txt", sep = ""))
  rsem <- read_table2(paste("../quant_r1/rsem/", dataset, "/1kg_NA12878_gencode100_", dataset , "_rsem.isoforms.results", sep = "")) %>%
    select(-effective_length, -expected_count, -TPM, -FPKM, -IsoPct)
  
  load(paste("../quant_r1/sim/", dataset, "/vg_r1/sim_1kg_NA12878_gencode100_", dataset , "_vg_r1.RData", sep = ""))
  
  truth_exp <- sim_exp %>%
    right_join(rsem, by = c("path" = "transcript_id")) %>%
    left_join(identical_seqs, by = c("path" = "Name1")) %>%
    replace_na(list(count = 0)) %>%
    group_by(path) %>% 
    mutate(n = n()) %>%
    mutate(count = count / n) %>%
    select(-n) %>%
    ungroup() %>%
    rename(name = Name2) %>%
    mutate(name = ifelse(is.na(name), paste(path, "_na", sep = ""), name)) %>%
    update_tpm() %>%
    group_by(name) %>%
    summarise(length = max(length), tpm_truth = sum(TPM), count_truth = sum(count))
  
  truth_exp <- truth_exp %>%
    ungroup() %>%
    mutate(count_truth = count_truth / sum(count_truth) * 1000000)
  
  if (read_type == "real") {
    
    truth_exp <- truth_exp %>% 
      mutate(length = 1) %>%
      mutate(tpm_truth = rep(c(1,2), length.out = nrow(truth_exp))) %>%
      mutate(count_truth = rep(c(1,2), length.out = nrow(truth_exp)))
  }
  
  print(sum(truth_exp$count_truth))
  print(sum(truth_exp$tpm_truth))
  
  # write.table(truth_exp, file = paste("truth/truth_exp_", read_type, "_", dataset, "_", ref_name, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
  
  files <- c(list.files(path = "methods", pattern = "rpvg.*.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "quant.sf.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "abundance.tsv.gz", full.names = T, recursive = T), list.files(path = "methods", pattern = "isoforms.results.gz", full.names = T, recursive = T))
  
  for (f in files) { 
    
    if (grepl("joint", f)) {
      
      next 
    }
    
    if (grepl("gibbs", f)) {
      
      next 
    }
    
    if (!grepl(dataset, f)) {
      
      next 
    }
    
    if (!grepl(read_type, f)) {
      
      next 
    }
    
    if (!grepl(paste("/", ref_name, sep = ""), f)) {
      
      next 
    }
    
    if (grepl("unidi", f) != grepl("unidi", ref_name)) {
      
      next 
    }
    
    print(f)
    
    if (grepl("quant.sf", f) | grepl("quant_boot.sf", f)) {
      
      exp_data <- parse_salmon(f)
      
    } else if (grepl("abundance.tsv", f)) {
      
      exp_data <- parse_kallisto(f)
      
    } else if (grepl("isoforms.results", f)) {
      
      exp_data <- parse_rsem(f)
      
    } else {
      
      exp_data <- parse_rpvg(f)
    }
    
    print(sum(exp_data$count_est))
    print(sum(exp_data$tpm_est))
    
    exp_data <- exp_data %>% 
      filter(name != "Unknown")
    
    if (read_type == "sim_vg") {
      
      exp_data <- exp_data %>%
        rename(count = count_est, TPM = tpm_est) %>%
        update_tpm() %>%
        rename(count_est = count, tpm_est = TPM)
    }
    
    exp_data <- exp_data %>% 
      ungroup() %>%
      mutate(count_est = count_est / sum(count_est) * 1000000)
    
    print(sum(exp_data$count_est))
    print(sum(exp_data$tpm_est))
    
    exp_data <- exp_data %>% 
      full_join(truth_exp, by = "name") %>%
      mutate(is_hap = ifelse(is.na(tpm_truth), F, T)) %>%
      replace_na(list(HaplotypeProbability = 0, tpm_est = 0, count_est = 0, Reads = exp_data$Reads[1], Method = exp_data$Method[1], Graph = exp_data$Graph[1], tpm_truth = 0, count_truth = 0)) %>%
      separate(name, c("transcript", "hap_id"), "_", extra = "drop", fill = "right") 
    
    exp_data_hap_prob <- exp_data %>%
      group_by(HaplotypeProbability, Reads, Method, Graph) %>%
      summarise(TP = sum((tpm_truth > 0) & (tpm_est > 0)),
                TN = sum((tpm_truth == 0) & (tpm_est == 0)),
                FP = sum((tpm_truth == 0) & (tpm_est > 0)),
                FN = sum((tpm_truth > 0) & (tpm_est == 0)),
                TP_tpm = sum((tpm_truth > 0) * tpm_est),
                FP_tpm = sum((tpm_truth == 0) * tpm_est))
    
    exp_data <- exp_data %>% 
      mutate(count_est = ifelse(HaplotypeProbability >= hap_prob_thres, count_est, 0)) %>%
      mutate(tpm_est = ifelse(HaplotypeProbability >= hap_prob_thres, tpm_est, 0))
    
    exp_data_stats_all <- exp_data %>%
      getStats() %>%
      add_column(Type = "All")
    
    exp_data_stats_exp <- exp_data %>%
      filter(tpm_est > 0) %>%
      getStats() %>%
      add_column(Type = "Expressed")
    
    exp_data_stats_hap <- exp_data %>%
      filter(is_hap) %>%
      getStats() %>%
      add_column(Type = "Haplotype")
    
    exp_data_stats_txp <- exp_data %>%
      group_by(transcript, Reads, Method, Graph) %>%
      summarise(count_est = sum(count_est), count_truth = sum(count_truth), tpm_est = sum(tpm_est), tpm_truth = sum(tpm_truth), is_hap = max(is_hap)) %>%
      ungroup() %>%
      getStats() %>%
      add_column(Type = "Transcript")
    
    exp_data_stats <- rbind(exp_data_stats_all, exp_data_stats_exp, exp_data_stats_hap, exp_data_stats_txp) %>%
      add_column(hap_prob_thres = hap_prob_thres) 
    
    exp_data_hap_exp <- exp_data %>%
      group_by(tpm_est, Reads, Method, Graph) %>%
      summarise(TP = sum((tpm_truth > 0) & (tpm_est > 0)),
                TN = sum((tpm_truth == 0) & (tpm_est == 0)),
                FP = sum((tpm_truth == 0) & (tpm_est > 0)),
                FN = sum((tpm_truth > 0) & (tpm_est == 0)),
                TP_tpm = sum((tpm_truth > 0) * tpm_est),
                FP_tpm = sum((tpm_truth == 0) * tpm_est)) %>%
      add_column(hap_prob_thres = hap_prob_thres) 
    
    exp_data_exp_is_hap <- exp_data %>%
      filter(is_hap) %>% 
      select(-transcript, -hap_id, -length.x, -HaplotypeProbability, -length.y, -is_hap) %>%
      add_column(hap_prob_thres = hap_prob_thres) 
    
    print(exp_data_stats)
    
    save(exp_data_hap_prob, exp_data_hap_exp, exp_data_exp_is_hap, exp_data_stats, file = paste("rdata/", exp_data_stats$Method[1], exp_data_stats$Reads[1], exp_data_stats$Graph[1], ".RData", sep = "", collapse = ""))
  }
}

dataset <- "SRR1153470"
sim_mean <- 277
sim_sd <- 43

for (ref_name in c("gencode100", "1kg_NA12878_gencode100", "1kg_EURnonCEU_af002_gencode100", "1kg_nonCEU_af001_gencode100", "1kg_all_af001_gencode100")) {
  
  parse_data(dataset, sim_mean, sim_sd, "sim_vg", ref_name)
}

for (ref_name in c("1kg_EURnonCEU_af002_gencode100", "1kg_EURnonCEU_af002_gencode100_unidi", "1kg_nonCEU_af001_gencode100", "1kg_nonCEU_af001_gencode100_unidi", "1kg_all_af001_gencode100", "1kg_all_af001_gencode100_unidi")) {
  
  parse_data(dataset, sim_mean, sim_sd, "real", ref_name)
}

dataset <- "ENCSR000AED_rep1"
sim_mean <- 216
sim_sd <- 24

for (ref_name in c("gencode100", "1kg_NA12878_gencode100", "1kg_EURnonCEU_af002_gencode100", "1kg_nonCEU_af001_gencode100", "1kg_all_af001_gencode100")) {
  
  parse_data(dataset, sim_mean, sim_sd, "sim_vg", ref_name)
}

for (ref_name in c("1kg_EURnonCEU_af002_gencode100", "1kg_EURnonCEU_af002_gencode100_unidi", "1kg_nonCEU_af001_gencode100", "1kg_nonCEU_af001_gencode100_unidi", "1kg_all_af001_gencode100", "1kg_all_af001_gencode100_unidi")) {
  
  parse_data(dataset, sim_mean, sim_sd, "real", ref_name)
}
