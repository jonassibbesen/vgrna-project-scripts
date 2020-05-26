
# plot_mapping_bias.R

rm(list=ls())

library("tidyverse")
library("gridExtra")

source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/plots/mapping/")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_tsv(filename, col_types = "icc")
  
  if (dir_split[2] == "vg") {
    
    data <- data %>%
      add_column(Method = dir_split[5]) %>%
      add_column(Graph = dir_split[6])
    
  } else {
    
    data <- data %>%
      add_column(Method = dir_split[2]) %>%
      add_column(Graph = dir_split[5])
  }
  
  data$Method = recode_factor(data$Method, "hisat2" = "HISAT2", "star" = "STAR", "map" = "vg map", "mpmap" = "vg mpmap", "mpmap_old" = "vg mpmap (old)")
  data$Graph = recode_factor(data$Graph, "gencode100" = "Linear", "1kg_nonCEU_af001_gencode100" = "1000g")
  
  return(data)
}

coverage_data_h1_0 <- map_dfr(list.files(pattern="*allele_cov_m0_h1.txt", full.names = T, recursive = T), parse_file)
coverage_data_h2_0 <- map_dfr(list.files(pattern="*allele_cov_m0_h2.txt", full.names = T, recursive = T), parse_file)

coverage_data_0 <- full_join(coverage_data_h1_0, coverage_data_h2_0, by = c("Position", "Method", "Graph"))

coverage_data_h1_20 <- map_dfr(list.files(pattern="*allele_cov_m20_h1.txt", full.names = T, recursive = T), parse_file)
coverage_data_h2_20 <- map_dfr(list.files(pattern="*allele_cov_m20_h2.txt", full.names = T, recursive = T), parse_file)

coverage_data_20 <- full_join(coverage_data_h1_20, coverage_data_h2_20, by = c("Position", "Method", "Graph"))

coverage_data_h1_40 <- map_dfr(list.files(pattern="*allele_cov_m40_h1.txt", full.names = T, recursive = T), parse_file)
coverage_data_h1_40 <- map_dfr(list.files(pattern="*allele_cov_m40_h2.txt", full.names = T, recursive = T), parse_file)

coverage_data_40 <- full_join(coverage_data_h1_40, coverage_data_h1_40, by = c("Position", "Method", "Graph"))

coverage_data <- rbind(coverage_data_0, coverage_data_20, coverage_data_40)
coverage_data$Position <- NULL

coverage_data <- coverage_data %>% 
  filter((AlleleType.x != AlleleType.y) & (AlleleType.x == "REF" | AlleleType.y == "REF")) %>% 
  rowwise() %>% 
  mutate(ref = ifelse(AlleleType.x == "REF", Count.x, Count.y)) %>%
  mutate(alt = ifelse(AlleleType.x != "REF", Count.x, Count.y)) %>%
  mutate(var = ifelse(AlleleType.x == "REF", AlleleType.y, AlleleType.x)) %>%
  filter(var != "COM") %>%
  mutate(binom_test = binom.test(c(Count.x, Count.y), p = 0.5, alternative = "two.sided")$p.value)

coverage_data$var <- factor(coverage_data$var, levels = c("SNV", "INS", "DEL"))
coverage_data$var = recode_factor(coverage_data$var, "SNV" = "SNVs", "INS" = "Insertions", "DEL" = "Deletions")

coverage_data %>%
  filter(Graph == "SNV") %>%
  ggplot(aes(x = binom_test)) +
  geom_histogram() +      
  facet_grid(Graph ~ Method) + 
  scale_y_continuous(limits=c(0, 500), oob = rescale_none) + 
  ggtitle("Petal and sepal length of iris")
  theme_bw() +
  theme(text = element_text(size=16))

coverage_data %>%
  filter(Graph == "Insertions") %>%
  ggplot(aes(x = binom_test)) +
  geom_histogram() +      
  facet_grid(Graph ~ Method) + 
  scale_y_continuous(limits=c(0, 500), oob = rescale_none)
  
coverage_data %>%
  filter(Graph == "Deletions") %>%
  ggplot(aes(x = binom_test)) +
  geom_histogram() +      
  facet_grid(Graph ~ Method) + 
  scale_y_continuous(limits=c(0, 500), oob = rescale_none)           
  
  ggplot(aes(sample = binom_test, col = Method)) +
  stat_qq()
  

pdf("mapping_bias.pdf", width = 10)

coverage_data %>% 
  filter((alt + ref) >= 10) %>%
  ggplot(aes(y = alt/(ref+alt), x = Method, fill = Graph)) +
  geom_violin() + 
  ylim(c(0,1)) +
  facet_grid(rows = vars(var)) + 
  ylab("Fraction mapped reads to alternative allele (MQ >= 40)") +
  labs(fill = "Methods", x = "") +
  theme_bw() +
  theme(text = element_text(size=16))

coverage_data %>% 
  filter((alt + ref) >= 10) %>%
  ggplot(aes(y = alt/(ref+alt), x = Method, fill = Graph)) +
  geom_violin() + 
  ylim(c(0.45,0.55)) +
  facet_grid(rows = vars(var)) + 
  ylab("Fraction mapped reads to alternative allele (MQ >= 40)") +
  labs(fill = "Methods", x = "") +
  theme_bw() +
  theme(text = element_text(size=16))

dev.off()


  


