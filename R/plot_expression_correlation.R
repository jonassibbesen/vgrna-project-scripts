
# plot_expression_correlation.R

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

parse_rpvg <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_names = T)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    rename(name = Name, tpm_est = TPM, count_est = ReadCount) %>%
    mutate(count_est = ifelse(HaplotypePosterior >= 0.95, count_est, 0)) %>%
    mutate(tpm_est = count_est / as.double(EffectiveLength)) %>%
    replace_na(list(tpm_est = 0)) %>%
    mutate(tpm_est = 10^6 * tpm_est / sum(tpm_est)) %>%
    select(-ClusterID, -Length, -HaplotypePosterior, -EffectiveLength, -ClusterRelativeExpression)
  
  return(data)
}

parse_kallisto <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_names = T)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    rename(name = target_id, tpm_est = tpm, count_est = est_counts) %>%
    select(-length, -eff_length)
  
  return(data)
}

parse_salmon <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_names = T)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9]) %>%
    rename(name = Name, tpm_est = TPM, count_est = NumReads) %>%
    select(-Length, -EffectiveLength)
  
  return(data)
}

#exp_data <- map_dfr(list.files(pattern = ".*rpvg_.*1kg_NA12878_gencode100_genes_sim.*.txt", full.names = T, recursive = T), parse_rpvg)
exp_data <- map_dfr(list.files(path = "./methods/kallisto/expression/polya_rna/sim_rsem/SRR1153470/kallisto/1kg_nonCEU_af001_gencode100/", pattern = ".*abundance.tsv", full.names = T, recursive = T), parse_kallisto)
#exp_data <- map_dfr(list.files(path = "./methods/salmon/expression/polya_rna/sim_rsem/SRR1153470/salmon/1kg_NA12878_gencode100/", pattern = ".*quant.sf", full.names = T, recursive = T), parse_salmon)

identical_seqs <- read_table2("graphs/1kg_NA12878_exons_gencode100_allpaths/1kg_nonCEU_af001_gencode100_hst_overlap.txt")

sim_exp_h1 <- read_table2("sim_1kg_NA12878_gencode100/SRR1153470/rsem/sim_1kg_NA12878_gencode100_SRR1153470_rsem_h1.sim.isoforms.results")
sim_exp_h2 <- read_table2("sim_1kg_NA12878_gencode100/SRR1153470/rsem/sim_1kg_NA12878_gencode100_SRR1153470_rsem_h2.sim.isoforms.results")

sim_exp <- bind_rows(sim_exp_h1, sim_exp_h2) %>%
  left_join(identical_seqs, by = c("transcript_id" = "Name1")) %>%
  rename(name = Name2) %>%
  mutate(name = ifelse(is.na(name), paste(transcript_id, "_na", sep = ""), name)) %>%
  mutate(TPM = count / effective_length) %>%
  replace_na(list(TPM = 0)) %>%
  mutate(TPM = 10^6 * TPM / sum(TPM)) %>%
  group_by(name) %>%
  summarise(tpm_sim = sum(TPM), count_sim = sum(count))

exp_data_sim <- exp_data %>% 
  full_join(sim_exp, by = "name") %>%
  mutate(is_hap = ifelse(is.na(tpm_sim), F, T)) %>%
  replace_na(list(HaplotypePosterior = 0, count_est = 0, tpm_est = 0, Reads = exp_data$Reads[1], Method = exp_data$Method[1], Graph = exp_data$Graph[1], tpm_sim = 0, count_sim = 0)) %>%
  select(-name, -count_est, -count_sim) 

print(nrow(exp_data_sim))

save(exp_data_sim, file = paste("rdata/", exp_data_sim$Method[1], exp_data_sim$Reads[1], exp_data_sim$Graph[1], ".RData", sep = "", collapse = ""))

