
# plot_coverage_correlation.R

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

coverage_data<- map_dfr(list.files(pattern=".*_real_.*_exon_cov_ENCSR706ANY.txt", full.names = T, recursive = T), parse_file)
pb_coverage <- read_table2("alignments/ENCSR706ANY/ENCSR706ANY_exon_cov_gc.txt", col_types = "iiciii")

coverage_data_mq30 <- coverage_data %>%
  mutate(Count = ifelse(MapQ < 0, 0, Count)) %>%
  mutate(ReadCoverage = Count * ReadCoverage) %>%
  mutate(BaseCoverage = Count * BaseCoverage) %>%
  group_by(AllelePosition, ExonSize, Reads, Method, Graph) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

pb_coverage <- pb_coverage %>%
  mutate(ReadCoverage = Count * ReadCoverage) %>%
  mutate(BaseCoverage = Count * BaseCoverage) %>%
  group_by(AllelePosition, ExonSize) %>%
  summarise(ReadCoverage = sum(ReadCoverage), BaseCoverage = sum(BaseCoverage))

coverage_data_pb_mq30 <- right_join(pb_coverage, coverage_data_mq30, by = c("AllelePosition", "ExonSize")) %>%
  mutate(BaseCoverage.x = BaseCoverage.x / ExonSize) %>%
  mutate(BaseCoverage.y = BaseCoverage.y / ExonSize)

coverage_data_pb_mq30 %>%
  group_by(Reads, Method, Graph) %>%
  summarise(n_base = sum(BaseCoverage.y * ExonSize) / 10^6, n_reads = sum(ReadCoverage.y), Pearson = cor(BaseCoverage.x, BaseCoverage.y, method = "pearson"), Spearman = cor(BaseCoverage.x, BaseCoverage.y, method = "spearman"))

coverage_data_pb_mq30 %>%
  filter(BaseCoverage.x > 0) %>%
  group_by(Reads, Method, Graph) %>%
  summarise(n_base = sum(BaseCoverage.y * ExonSize) / 10^6, n_reads = sum(ReadCoverage.y), Pearson = cor(BaseCoverage.x, BaseCoverage.y, method = "pearson"), Spearman = cor(BaseCoverage.x, BaseCoverage.y, method = "spearman"))


pdf("scatter.pdf")
coverage_data_pb_mq30 %>%
  filter(Graph == "1kg_nonCEU_af001_gencode100") %>%
  filter(BaseCoverage.x > 0 | BaseCoverage.y > 0) %>%
  group_by(Reads, Method, Graph) %>%
  mutate(BaseCoverage.x = BaseCoverage.x / sum(BaseCoverage.x) * 10^6) %>%
  mutate(BaseCoverage.y = BaseCoverage.y / sum(BaseCoverage.y) * 10^6) %>%
  ggplot(aes(x = log2(BaseCoverage.x + 1), y = log2(BaseCoverage.y + 1) - log2(BaseCoverage.x + 1))) +
  geom_abline(intercept = 0, color = "red") +
  geom_point(alpha = 0.1) +
  facet_grid(. ~ Method) +
  scale_fill_brewer(palette="Dark2") +
  coord_fixed() +
  labs(fill = "Methods", x = "") +
  xlab("log10(1 + simulated TPM)") +
  ylab("log10(1 + estimated TPM)") +
  theme_bw() +
  theme(strip.background = element_blank(), panel.border = element_rect(colour = "black")) +
  theme(legend.position = "none") +
  theme(text = element_text(size=14))
dev.off()




