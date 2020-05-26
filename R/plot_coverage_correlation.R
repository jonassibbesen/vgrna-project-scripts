
# plot_overlap_benchmark.R
# Plot mapping coverage benchmark 

rm(list=ls())

library("tidyverse")
library("gridExtra")

# args <- commandArgs()
# script_dir <- dirname(sub("--file=", "", args[4]))
# print(script_dir)
# 
# print(args)
# system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
# system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))
# 
# data_dir <- read.csv(args[6], sep = " ", header = F)

setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/plots/debug/mapping3/1kg_NA12878_exons_gencode/real_SRR1153470/")

parse_file <- function(filename, pb_coverage) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_tsv(filename, col_names = F, col_types = "ciii")
  data <- data %>%
    right_join(pb_coverage, by = c("X1", "X2", "X3")) %>%
    add_column(Method = dir_split[2]) 
  
  print(nrow(data))
  return(data)
}

pb_coverage <- read_tsv("pacbio_coverage_mq40.txt", col_names = F, col_types = "ciii")

coverage_data <- map_dfr(list.files(pattern="*_5M_coverage_mq40.txt", full.names = T, recursive = T), parse_file, pb_coverage)

coverage_data <- coverage_data %>%
  filter(Method != "mpmap")

coverage_data$Method = recode_factor(coverage_data$Method, "hisat2" = "Hisat2", "map" = "vg map", "mpmap_mapq" = "vg mpmap")

#coverage_data$Method = recode_factor(coverage_data$Method, "hisat2" = "Hisat2", "map" = "vg map", "mpmap" = "vg mpmap", "mpmap_mapq" = "vg mpmap (new mapq)")

coverage_data_cor <- coverage_data %>%
  group_by(Method) %>%
  mutate(pearson = cor(X4.x, X4.y, method = "pearson")) %>%
  mutate(spearman = cor(X4.x, X4.y, method = "spearman")) %>%
  summarise(n = n(), Pearson = max(pearson), Spearman = max(spearman), sum.x = sum(X4.x / 10^9), sum.y = sum(X4.y / 10^9))

pdf("coverage_cor.pdf", width = 4)
ggplot(coverage_data_cor, aes(y = Pearson, x = Method, fill = Method)) +
  geom_bar(position = position_dodge2(preserve = "single"), stat="identity") +
  labs(fill = "Methods", x = "") +
  ylab("PacBio exon coverage Pearson correlation") +
  theme_bw() +
  theme(legend.position="none") +
  theme(text = element_text(size=20))
dev.off()

pdf("coverage_base.pdf", width = 4)
ggplot(coverage_data_cor, aes(y = sum.x, x = Method, fill = Method)) +
  geom_bar(position = position_dodge2(preserve = "single"), stat="identity") +
  labs(fill = "Methods", x = "") +
  ylab("Number of aligned exon bases (10^9)") +
  theme_bw() +
  theme(legend.position="none") +
  theme(text = element_text(size=20))
dev.off()

pdf("scat.pdf")
coverage_data %>%
  filter(X1 == "1") %>%
  filter(Method == "vg mpmap") %>%
  ggplot(aes(x = log10(X4.x + 1), y = log10(X4.y + 1))) +
  geom_point(alpha = 0.1) +
  facet_grid(. ~ Method) +
  scale_fill_brewer(palette="Dark2") +
  coord_fixed() +
  labs(fill = "Methods", x = "") +
  xlab("log10(1 + vg mpmap exon base coverage)") +
  ylab("log10(1 + PacBio exon base coverage)") +
  theme_bw() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black")) +
  theme(legend.position = "none") +
  theme(text = element_text(size=16))
dev.off()