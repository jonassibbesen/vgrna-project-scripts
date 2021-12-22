
# plot_allele_specific_expression.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("scales")
library("wesanderson")
library("truncnorm")

# source("./utils.R")

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/ase_r1/")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

########

p_threshold <- 10^-4

#dataset <- "SRR1153470"
#dataset_wasp <- "SRR"

dataset <- "ENCSR000AED_rep1"
dataset_wasp <- "ENC"

print(dataset)
print(dataset_wasp)

parse_het_vars <- function(filename) {
  
  data <- read_table2(filename, col_names = F, col_types = "ciccc")
  return(data)
}

parse_ase_data <- function(filename, het_vars) {
  
  data <- read_table2(filename, col_names = T, col_types = "ciccciiidddd") 
  return(data)
}

### 

het_vars <- map_dfr(list.files(path = "het_vars", pattern = "1kg_NA12878_exons.*_het.txt", full.names = T, recursive = T), parse_het_vars) %>%
  separate(X5, c("GT1", "GT2"), "\\|") 

colnames(het_vars) <- c("Chrom", "Position", "Ref", "Alt", "GT1", "GT2")

print(het_vars %>% filter(GT1 == GT2))

het_vars <- het_vars %>% 
  rowwise() %>%
  mutate(Allele1 = str_split(paste(Ref, Alt, sep = ","), ",")[[1]][as.integer(GT1) + 1]) %>%
  mutate(Allele2 = str_split(paste(Ref, Alt, sep = ","), ",")[[1]][as.integer(GT2) + 1]) 

het_vars <-het_vars %>%
  select(-Ref, -Alt, -GT1, -GT2)

ase_data_sim <- map_dfr(list.files(path = paste("allele_expression/", dataset, "/sim_vg_r1/1kg_NA12878_gencode100", sep = ""), pattern = "allele_exp_sim_vg.*.txt", full.names = T, recursive = T), parse_ase_data) %>%
  select(-AlleleNum, -AlleleLength, -HomopolymerLength, -NumTandemRepeats, -Probability, -TranscriptReadCount, -TPM) %>%
  rename(base_count_sim = BaseReadCount)

ase_base_count_sim <- het_vars %>% left_join(ase_data_sim, by = c("Chrom", "Position", "Allele1" = "AlleleSeq"))
ase_base_count_sim <- ase_base_count_sim %>% left_join(ase_data_sim, by = c("Chrom", "Position", "Allele2" = "AlleleSeq"), suffix = c(".1", ".2")) 

ase_base_count_sim %>% filter(is.na(base_count_sim.1))
ase_base_count_sim %>% filter(is.na(base_count_sim.2))

ase_base_count_sim <- ase_base_count_sim %>%
  rowwise() %>%
  mutate(binom_test_sim = binom.test(x = c(round(base_count_sim.1), round(base_count_sim.2)), alternative = "two.sided")$p.value)

### 

ase_data_wasp <- read_table2(paste("wasp/as_counts_", dataset_wasp, "_norm.txt.gz", sep = ""), col_names = F, col_types = "ciccciii") %>%
  separate(X5, c("GT1", "GT2"), "\\|") %>%
  filter(GT1 != GT2) %>%
  select(-X8) %>%
  rowwise() %>%
  mutate(X1 = substr(X1, 4, nchar(X1)))

colnames(ase_data_wasp) <- c("Chrom", "Position", "Ref", "Alt", "GT1", "GT2", "base_count_est.1", "base_count_est.2")

ase_data_wasp <- ase_data_wasp %>%
  rowwise() %>%
  mutate(Allele1 = str_split(paste(Ref, Alt, sep = ","), ",")[[1]][as.integer(GT1) + 1]) %>%
  mutate(Allele2 = str_split(paste(Ref, Alt, sep = ","), ",")[[1]][as.integer(GT2) + 1])

ase_data_wasp <-ase_data_wasp %>%
  select(-Ref, -Alt, -GT1, -GT2)

ase_base_count_sim_wasp <- ase_base_count_sim %>% left_join(ase_data_wasp, by = c("Chrom", "Position", "Allele1", "Allele2"))

ase_base_count_sim_wasp <- ase_base_count_sim_wasp %>%
  replace_na(list(base_count_est.1 = 0, base_count_est.2 = 0)) %>%
  mutate(Method = "WASP (STAR)") %>%
  mutate(Variants = "1kg_NA12878") %>%
  select(Method, Variants, AlleleType.1, AlleleType.2, base_count_sim.1, base_count_sim.2, binom_test_sim, base_count_est.1, base_count_est.2)

ase_base_count_sim_est <- ase_base_count_sim_wasp

### 

ase_data_rpvg_na <- map_dfr(list.files(path = paste("allele_expression/", dataset, "/inference_r1/", sep = ""), pattern = paste("allele_exp_rpvg_mpmap_1kg_NA12878", ".*.txt", sep = ""), full.names = T, recursive = T), parse_ase_data) %>%
  select(-AlleleNum, -AlleleType, -AlleleLength, -HomopolymerLength, -NumTandemRepeats, -Probability, -TranscriptReadCount, -TPM) %>%
  rename(base_count_est = BaseReadCount) 

ase_base_count_sim_rpvg_na <- ase_base_count_sim %>% left_join(ase_data_rpvg_na, by = c("Chrom", "Position", "Allele1" = "AlleleSeq"))
ase_base_count_sim_rpvg_na <- ase_base_count_sim_rpvg_na %>% left_join(ase_data_rpvg_na, by = c("Chrom", "Position", "Allele2" = "AlleleSeq"), suffix = c(".1", ".2"))

ase_base_count_sim_rpvg_na <- ase_base_count_sim_rpvg_na %>%
  replace_na(list(base_count_est.1 = 0, base_count_est.2 = 0)) %>%
  mutate(Method = "mpmap-rpvg") %>%
  mutate(Variants = "1kg_NA12878") %>%
  select(Method, Variants, AlleleType.1, AlleleType.2, base_count_sim.1, base_count_sim.2, binom_test_sim, base_count_est.1, base_count_est.2)

ase_base_count_sim_est <- rbind(ase_base_count_sim_est, ase_base_count_sim_rpvg_na)

### 
# 
# for (var in c("1kg_EURnonCEU")) {
#   
#   ase_data_rpvg <- map_dfr(list.files(path = "allele_expression/SRR1153470/inference_r1/", pattern = paste("allele_exp_rpvg_mpmap_", var, ".*.txt", sep = ""), full.names = T, recursive = T), parse_ase_data) %>%
#     select(-AlleleNum, -AlleleType, -AlleleLength, -HomopolymerLength, -NumTandemRepeats, -Probability, -TranscriptReadCount, -TPM) %>%
#     rename(base_count_est = BaseReadCount) %>%
#     group_by(Chrom, Position) %>%
#     mutate(n_var_exp = sum(base_count_est > 0)) %>%
#     filter(n_var_exp >= 2) 
#     
#   ase_base_count_sim_rpvg <- ase_base_count_sim %>% left_join(ase_data_rpvg, by = c("Chrom", "Position", "Allele1" = "AlleleSeq"))
#   ase_base_count_sim_rpvg <- ase_base_count_sim_rpvg %>% left_join(ase_data_rpvg, by = c("Chrom", "Position", "Allele2" = "AlleleSeq"), suffix = c(".1", ".2")) %>%
#     replace_na(list(base_count_est.1 = 0, base_count_est.2 = 0))
#   
#   ase_base_count_sim_rpvg <- ase_base_count_sim_rpvg %>%
#     mutate(Method = "mpmap-rpvg") %>%
#     mutate(Variants = var) %>%
#     select(Method, Variants, AlleleType.1, AlleleType.2, base_count_sim.1, base_count_sim.2, binom_test_sim, base_count_est.1, base_count_est.2)
#   
#   ase_base_count_sim_est <- rbind(ase_base_count_sim_est, ase_base_count_sim_rpvg)
# }

### 

ase_base_count_sim_est <- ase_base_count_sim_est %>%
  add_column(Count = 1) %>%
  group_by(Method, Variants, AlleleType.1, AlleleType.2, base_count_sim.1, base_count_sim.2, binom_test_sim, base_count_est.1, base_count_est.2) %>%
  summarise(Count = sum(Count))

ase_base_count_sim_est <- ase_base_count_sim_est %>%
  rowwise() %>%
  mutate(binom_test_est = binom.test(x = c(round(base_count_est.1), round(base_count_est.2)), alternative = "two.sided")$p.value)

ase_base_count_sim_est <- ase_base_count_sim_est %>%
  rowwise() %>%
  mutate(VarType = ifelse(AlleleType.1 == "Ref", AlleleType.2, AlleleType.1)) %>%
  mutate(VarType = ifelse(AlleleType.1 != "Ref" & AlleleType.2 != "Ref" & AlleleType.1 == AlleleType.2, AlleleType.1, VarType))  %>%
  mutate(VarType = ifelse(AlleleType.1 != "Ref" & AlleleType.2 != "Ref" & AlleleType.1 != AlleleType.2, "Multi", VarType))


ase_base_count_sim_est %>%
  filter(base_count_sim.1 > 0 | base_count_sim.2 > 0) %>%
  filter(VarType == "SNV" | VarType == "Del" | VarType == "Ins") %>%
  group_by(Method, Variants, VarType) %>%
  summarise(TPR = sum(Count * (binom_test_sim <= p_threshold & binom_test_est <= p_threshold)) / sum(Count * (binom_test_sim <= p_threshold)), 
            FPR = sum(Count * (binom_test_est <= p_threshold & binom_test_sim > p_threshold)) / sum(Count),
            n_exp_est = sum(Count * (base_count_est.1 + base_count_est.2 > 0)), 
            n_ase_est = sum(Count * (binom_test_est <= p_threshold)), 
            n_ase_sim = sum(Count * (binom_test_sim <= p_threshold)), 
            n_total = sum(Count)) %>%
  print(n=100)

ase_base_count_sim_est_pval_roc <- ase_base_count_sim_est %>%
  filter(base_count_sim.1 > 0 | base_count_sim.2 > 0) %>%
  filter(VarType == "SNV" | VarType == "Del" | VarType == "Ins") %>%
  mutate(var_exp_est = base_count_est.1 + base_count_est.2) %>%
  mutate(Pn = (Count *(binom_test_sim <= p_threshold))) %>%
  mutate(TP = (Count * (binom_test_sim <= p_threshold & binom_test_est <= p_threshold))) %>%
  mutate(FP = (Count * (binom_test_sim > p_threshold & binom_test_est <= p_threshold))) %>%
  group_by(Method, Variants, VarType, var_exp_est) %>%
  summarize(Pn = sum(Pn), TP = sum(TP), FP = sum(FP), Count = sum(Count)) %>%
  arrange(desc(var_exp_est), .by_group = T) %>%
  mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
  group_by(Method, Variants, VarType) %>%
  mutate(TPR = TPcs / sum(Pn), Pn = sum(Pn), FPR = FPcs / sum(Count), Count = sum(Count))

print(ase_base_count_sim_est_pval_roc %>% select(-TP, -FP) %>% tail())

ase_base_count_sim_est_pval_roc$Method <- factor(ase_base_count_sim_est_pval_roc$Method, levels = c("WASP (STAR)", "mpmap-rpvg"))

ase_base_count_sim_est_pval_roc$VarType = recode_factor(ase_base_count_sim_est_pval_roc$VarType, 
                                                     "SNV" = "SNV", 
                                                     "Ins" = "Insertion", 
                                                     "Del" = "Deletion")

ase_base_count_sim_est_pval_roc$VarType <- factor(ase_base_count_sim_est_pval_roc$VarType, levels = c("SNV", "Insertion", "Deletion"))

ase_base_count_sim_est_pval_roc$Variants = recode_factor(ase_base_count_sim_est_pval_roc$Variants,
                                                         "1kg_NA12878" = "Sample-specific (NA12878)",
                                                         "1kg_EURnonCEU" = "Europe (excl. CEU)",
                                                         "1kg_nonCEU" = "Whole (excl. CEU)",
                                                         "1kg_all" = "Whole")

ase_base_count_sim_est_pval_roc$Variants <- factor(ase_base_count_sim_est_pval_roc$Variants, levels = c("Sample-specific (NA12878)", "Europe (excl. CEU)", "Whole (excl. CEU)", "Whole"))

ase_base_count_sim_est_pval_roc$FacetRow <- "Simulated reads (vg)"

wes_cols <- c(wes_palette("GrandBudapest1")[1], wes_palette("Chevalier1")[1])

pdf(paste("plots/sim_r1_ase_roc_", dataset, ".pdf", sep = ""), height = 5, width = 9, pointsize = 12)
ase_base_count_sim_est_pval_roc %>%
  ggplot(aes(y = TPR, x = FPR, color = Method, linetype = Variants, shape = Variants)) +
  geom_line(size = 0.75) +
  scale_color_manual(values = wes_cols) +
  facet_grid(FacetRow ~ VarType) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ylim(c(0,1)) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(0.5, "cm")) +
  theme(text = element_text(size = 10))
dev.off()

