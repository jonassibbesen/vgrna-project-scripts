
# plot_region_overlap.R

rm(list=ls())

library("tidyverse")
library("gridExtra")
library("wesanderson")
library("scales")

# source("./utils.R")

# printHeader()

# data_dir <- read.csv(args[6], sep = " ", header = F)
# setwd(data_dir)

source("/Users/jonas/Documents/postdoc/sc/code/vgrna-project-scripts/R/utils.R")
setwd("/Users/jonas/Documents/postdoc/sc/projects/vgrna/figures/mapping/")

parse_file <- function(filename) {
  
  dir_split <- strsplit(dirname(filename), "/")[[1]]
  
  data <- read_table2(filename)
  data <- data %>%
    add_column(Reads = dir_split[7]) %>%
    add_column(Method = dir_split[8]) %>%
    add_column(Graph = dir_split[9])
  
  if (grepl("ovl_ENCSR706ANY_mq0", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "ENCSR706ANY")
    
  } else if (grepl("ovl_ENCSR706ANY_mq30", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "ENCSR706ANY  (MapQ >= 30)")

  } else if (grepl("ovl_gc", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "GENCODE")

  } else {
    
    stopifnot(FALSE)
  }
  
  return(data)
}

overlap_data_raw <- map_dfr(list.files(path = "./methods", pattern=".*_exon_ovl_.*.txt", full.names = T, recursive = T), parse_file)

overlap_data_o05 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.05) %>%
  add_column(Threshold = "Overlap >= 5%") 

overlap_data_o75 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.75) %>%
  add_column(Threshold = "Overlap >= 75%") 

overlap_data <- bind_rows(overlap_data_o05, overlap_data_o75)


########


overlap_data_polya <- overlap_data %>%
  filter(Reads == data_set3) 

overlap_data_polya$Method <- recode_factor(overlap_data_polya$Method, 
                                           "hisat2" = "HISAT2",
                                           "star" = "STAR",
                                           "map_fast" = "vg map", 
                                           "mpmap" = "vg mpmap")

overlap_data_polya$Reads <- recode_factor(overlap_data_polya$Reads, 
                                          "SRR1153470_uni" = "Training set", 
                                          "ENCSR000AED_rep1_uni" = "Test set")

overlap_data_polya$Graph = recode_factor(overlap_data_polya$Graph, 
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)",
                                         "1kg_nonCEU_af001_gencode100_gtex10s2r8e1g" = "1000g (GTEx)")

overlap_data_polya_main <- overlap_data_polya

overlap_data_polya_main$FacetRow <- recode_factor(overlap_data_polya_main$Threshold, 
                                       "Overlap >= 75%" = "Overlap >= 75%", 
                                       "Overlap >= 5%" = "Overlap >= 5%")

overlap_data_polya_main$FacetCol <- recode_factor(overlap_data_polya_main$Regions, 
                                       "ENCSR706ANY" = "ENCSR706ANY", 
                                       "ENCSR706ANY  (MapQ >= 30)" = "ENCSR706ANY  (MapQ >= 30)", 
                                       "GENCODE" = "GENCODE")

plotOverlapBenchmarkMapQ(overlap_data_polya_main, wes_cols, "plots/polya_rna/real_overlap_polya_main")


overlap_data_polya_bar <- overlap_data_polya %>%
  filter(Regions == "GENCODE") %>%
  filter(Threshold == "Overlap >= 5%") %>%
  mutate(FacetRow = "") %>%
  mutate(FacetCol = "") %>%
  mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
  mutate(MapQ0 = Count * (MapQ >= 0)) %>% 
  mutate(MapQ1 = Count * (MapQ >= 1)) %>% 
  group_by(Method, Graph) %>%
  summarise(count = sum(Count), MapQ0 = sum(MapQ0), MapQ1 = sum(MapQ1)) %>%
  mutate(MapQ0_frac = MapQ0 / count) %>%
  mutate(MapQ1_frac = MapQ1 / count) %>%
  gather("MapQ0_frac", "MapQ1_frac", key = "Filter", value = "Frac")

overlap_data_polya_bar <- overlap_data_polya_bar %>%
 ungroup() %>%
  add_row(Method = "STAR", Graph = "1000g (no-CEU)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
  add_row(Method = "STAR", Graph = "1000g (no-CEU)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0)

overlap_data_polya_bar$Filter <- recode_factor(overlap_data_polya_bar$Filter, 
                                                "MapQ0_frac" = "All",
                                                "MapQ1_frac" = "MapQ > 0")

overlap_data_polya_bar$Graph <- recode_factor(overlap_data_polya_bar$Graph, 
                                            "Spliced reference" = "Spliced\nreference",
                                            "1000g (no-CEU)" = "1000g\n(no-CEU)")

pdf("plots/polya_rna/real_overlap_polya_bar.pdf", height = 4, width = 4, pointsize = 12)
ggplot() +
  geom_bar(data = overlap_data_polya_bar[overlap_data_polya_bar$Filter == "MapQ > 0",], aes(Graph, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = overlap_data_polya_bar[overlap_data_polya_bar$Filter == "All",], aes(Graph, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
  scale_fill_manual(values = wes_cols) +
  scale_alpha_manual(name = "Filter", values = c(0.5, 1), labels = c("Unfiltered", "MapQ > 0"), drop = F) + 
  scale_y_continuous(limits=c(0.85, 1), oob = rescale_none) +
  xlab("") +
  ylab("Mapping rate") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=14))
dev.off()


########
