
rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping/")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename, col_types = "iiciii")
  data <- data %>%
    add_column(Reads = dir_split[6]) %>%
    add_column(Method = dir_split[7]) %>%
    add_column(Graph = dir_split[8])
  
  return(data)
}

coverage_data<- map_dfr(list.files(pattern=".*_real_.*_exon_cov.txt", full.names = T, recursive = T), parse_file)
pb_coverage <- read_table2("ENCSR706ANY/ENCSR706ANY_exon_cov.txt", col_types = "iiciii")

coverage_data_mq30 <- coverage_data %>%
  filter(MapQ > 0) %>%
  group_by(AllelePosition, ExonSize, Reads, Method, Graph) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

pb_coverage_mq30 <- pb_coverage %>%
  filter(MapQ > 0) %>%
  group_by(AllelePosition, ExonSize) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

coverage_data_pb_mq30 <- full_join(pb_coverage_mq30, coverage_data_mq30, by = c("AllelePosition", "ExonSize")) %>%
  replace_na(list(Count.x = 1, Count.y = 1)) %>%
  replace_na(list(ReadCoverage.x = 0, ReadCoverage.y = 0)) %>%
  replace_na(list(BaseCoverage.x = 0, BaseCoverage.y = 0)) %>%
  mutate(BaseCoverage.x = BaseCoverage.x / ExonSize) %>%
  mutate(BaseCoverage.y = BaseCoverage.y / ExonSize)

coverage_data_pb_mq30_cor <- coverage_data_pb_mq30 %>%
  group_by(Reads, Method, Graph) %>%
  mutate(BaseCoverage.x_norm = BaseCoverage.x / sum(BaseCoverage.x) * 10^6) %>%
  mutate(BaseCoverage.y_norm = BaseCoverage.y / sum(BaseCoverage.y) * 10^6) %>%
  mutate(ard = abs(BaseCoverage.x_norm - BaseCoverage.y_norm) / (BaseCoverage.x_norm + BaseCoverage.y_norm)) %>%
  replace_na(list(ard = 0)) %>%
  summarise(n_base = sum(BaseCoverage.x * ExonSize) / 10^6, n_reads = sum(ReadCoverage.x), Pearson = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "pearson"), Spearman = cor(BaseCoverage.x_norm, BaseCoverage.y_norm, method = "spearman"), ard_mean = mean(ard), ard_median = median(ard))




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