
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

  } else if (grepl("ovl_mb", basename(filename))) {
    
    data <- data %>%
      add_column(Regions = "miRBase")
        
  } else {
    
    stopifnot(FALSE)
  }
  
  return(data)
}

overlap_data_raw <- map_dfr(list.files(pattern=".*_exon_ovl_.*.txt", full.names = T, recursive = T), parse_file)

overlap_data_o10 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.1) %>%
  add_column(Threshold = "Overlap >= 10%") 

overlap_data_o90 <- overlap_data_raw %>%
  mutate(Correct = Overlap >= 0.9) %>%
  add_column(Threshold = "Overlap >= 90%") 

overlap_data <- bind_rows(overlap_data_o10, overlap_data_o90) 


########


overlap_data_polya <- overlap_data %>%
  filter(Reads == data_set3) 

overlap_data_polya$Method <- recode_factor(overlap_data_polya$Method, 
                                           "hisat2" = "HISAT2",
                                           "star" = "STAR",
                                           "map" = "vg map (def)", 
                                           "map_fast" = "vg map", 
                                           "mpmap" = "vg mpmap")

overlap_data_polya <- overlap_data_polya %>%
  filter(Method != "vg map (def)")

overlap_data_polya$Reads <- recode_factor(overlap_data_polya$Reads, 
                                          "SRR1153470_uni" = "Training set", 
                                          "ENCSR000AED_rep1_uni" = "Test set")

overlap_data_polya$Graph = recode_factor(overlap_data_polya$Graph, 
                                         "gencode100" = "Spliced reference",
                                         "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)",
                                         "1kg_nonCEU_af001_gencode100_gtex10s2r8e1g" = "1000g (GTEx)")

overlap_data_polya_main <- overlap_data_polya

overlap_data_polya_main$FacetRow <- recode_factor(overlap_data_polya_main$Threshold, 
                                       "Overlap >= 90%" = "Overlap >= 90%", 
                                       "Overlap >= 10%" = "Overlap >= 10%")

overlap_data_polya_main$FacetCol <- recode_factor(overlap_data_polya_main$Regions, 
                                       "ENCSR706ANY" = "ENCSR706ANY", 
                                       "ENCSR706ANY  (MapQ >= 30)" = "ENCSR706ANY  (MapQ >= 30)", 
                                       "GENCODE" = "GENCODE")

plotOverlapBenchmark(overlap_data_polya_main, wes_cols, "plots/polya_rna/real_overlap_polya_main", F)


overlap_data_polya_mapq <- overlap_data_polya %>%
  filter(Regions == "GENCODE") %>%
  filter(Threshold == "Overlap >= 10%") %>%
  mutate(FacetRow = "") %>%
  mutate(FacetCol = "") 

plotMapQBenchmark(overlap_data_polya_mapq, wes_cols, "plots/polya_rna/real_overlap_polya_mapq")


overlap_data_polya_mapq_bar <- overlap_data_polya_mapq %>% 
  mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
  mutate(MapQ0 = Count * (MapQ >= 0)) %>% 
  mutate(MapQ1 = Count * (MapQ >= 1)) %>% 
  group_by(Method, Graph) %>%
  summarise(count = sum(Count), MapQ0 = sum(MapQ0), MapQ1 = sum(MapQ1)) %>%
  mutate(MapQ0_frac = MapQ0 / count) %>%
  mutate(MapQ1_frac = MapQ1 / count) %>%
  gather("MapQ0_frac", "MapQ1_frac", key = "Filter", value = "Frac")

overlap_data_polya_mapq_bar <- overlap_data_polya_mapq_bar %>%
 ungroup() %>%
  add_row(Method = "STAR", Graph = "1000g (no-CEU)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
  add_row(Method = "STAR", Graph = "1000g (no-CEU)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0) %>%
  add_row(Method = "HISAT2", Graph = "1000g (GTEx)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
  add_row(Method = "HISAT2", Graph = "1000g (GTEx)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0) %>%
  add_row(Method = "STAR", Graph = "1000g (GTEx)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
  add_row(Method = "STAR", Graph = "1000g (GTEx)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0) 
  
overlap_data_polya_mapq_bar$Filter = recode_factor(overlap_data_polya_mapq_bar$Filter, 
                                                "MapQ0_frac" = "All",
                                                "MapQ1_frac" = "MapQ > 0")

overlap_data_polya_mapq_bar$Graph = recode_factor(overlap_data_polya_mapq_bar$Graph, 
                                            "Spliced reference" = "Spliced\nreference",
                                            "1000g (no-CEU)" = "1000g\n(no-CEU)",
                                            "1000g (GTEx)" = "1000g\n(GTEx)")

pdf("plots/polya_rna/real_overlap_polya_mapq_bar.pdf", height = 4, width = 4.5, pointsize = 12)
ggplot() +
  geom_bar(data = overlap_data_polya_mapq_bar[overlap_data_polya_mapq_bar$Filter == "MapQ > 0",], aes(Graph, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge()) +
  geom_bar(data = overlap_data_polya_mapq_bar[overlap_data_polya_mapq_bar$Filter == "All",], aes(Graph, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
  scale_fill_manual(values = wes_cols) +
  scale_alpha_manual(name = "Filter", values = c(0.5, 1), labels = c("Unfiltered", "MapQ > 0"), drop = F) + 
  scale_y_continuous(limits=c(0.85, 1), oob = rescale_none) +
  xlab("") +
  ylab("Mapping rate") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=12))
dev.off()


########

# 
# overlap_data_mir <- overlap_data %>%
#   filter(Reads == "ENCSR958UOC_rep1") 
# 
# overlap_data_mir$Method <- recode_factor(overlap_data_mir$Method, 
#                                          "bowtie2" = "Bowtie2",
#                                          "bowtie2_vs_end" = "Bowtie2 (vs end)",
#                                          "bowtie2_vs_local" = "Bowtie2 (vs local)",
#                                          "hisat2_nosplice" = "HISAT2", 
#                                          "star_nosplice" = "STAR", 
#                                          "star_encode" = "STAR (encode)",
#                                          "map" = "vg map (def)", 
#                                          "map_fast" = "vg map", 
#                                          "mpmap" = "vg mpmap (multi)", 
#                                          "mpmap_nomulti" = "vg mpmap")
# 
# overlap_data_mir <- overlap_data_mir %>%
#   filter(Method != "Bowtie2 (vs end)") %>%
#   filter(Method != "Bowtie2 (vs local)") %>%
#   filter(Method != "STAR") %>% 
#   filter(Method != "vg map (def)") %>%
#   filter(Method != "vg mpmap (multi)")
# 
# overlap_data_mir$Reads <- recode_factor(overlap_data_mir$Reads, 
#                                         "ERR187607" = "Training set (old)", 
#                                         "ENCSR958UOC_rep1" = "Training set (rep1)", 
#                                         "ENCSR958UOC_rep2" = "Training set (rep2)")
# 
# overlap_data_mir$Graph = recode_factor(overlap_data_mir$Graph, 
#                                        "linear" = "Linear",
#                                        "gencode100" = "Linear",
#                                        "1kg_nonCEU_af001_gencode100" = "1000g (no-CEU)")
# 
#   
# 
# wes_cols_mir <- wes_cols[c(6, seq(1, 5))]
# 
# overlap_data_mir_main <- overlap_data_mir
# 
# overlap_data_mir_main$FacetCol <- recode_factor(overlap_data_mir_main$Threshold, 
#                                                   "Overlap >= 90%" = "Overlap >= 90%", 
#                                                   "Overlap >= 10%" = "Overlap >= 10%")
# 
# overlap_data_mir_main$FacetRow <- ""
# 
# plotOverlapBenchmark(overlap_data_mir_main, wes_cols_mir, "plots/micro_rna/real_overlap_mir_main", F)
# 
# 
# overlap_data_mir_mapq <- overlap_data_mir %>%
#   filter(Regions == "miRBase") %>%
#   filter(Threshold == "Overlap >= 10%") %>%
#   mutate(FacetRow = "") %>%
#   mutate(FacetCol = "") 
# 
# plotMapQBenchmark(overlap_data_mir_mapq, wes_cols_mir, "plots/micro_rna/real_overlap_mir_mapq")
# 
# 
# overlap_data_mir_mapq_bar <- overlap_data_mir_mapq %>% 
#   mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
#   mutate(MapQ0 = Count * (MapQ >= 0)) %>% 
#   mutate(MapQ1 = Count * (MapQ >= 1)) %>% 
#   group_by(Method, Graph, Reads) %>%
#   summarise(count = sum(Count), MapQ0 = sum(MapQ0), MapQ1 = sum(MapQ1)) %>%
#   mutate(MapQ0_frac = MapQ0 / count) %>%
#   mutate(MapQ1_frac = MapQ1 / count) %>%
#   gather("MapQ0_frac", "MapQ1_frac", key = "Filter", value = "Frac") 
# 
# overlap_data_mir_mapq_bar <- overlap_data_mir_mapq_bar %>%
#   ungroup() %>%
#   add_row(Method = "Bowtie2", Graph = "1000g (no-CEU)", Reads = "Training set (rep1)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
#   add_row(Method = "Bowtie2", Graph = "1000g (no-CEU)", Reads = "Training set (rep1)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0) %>%
#   add_row(Method = "Bowtie2", Graph = "1000g (no-CEU)", Reads = "Training set (rep2)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
#   add_row(Method = "Bowtie2", Graph = "1000g (no-CEU)", Reads = "Training set (rep2)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0) %>%
#   add_row(Method = "STAR (encode)", Graph = "1000g (no-CEU)", Reads = "Training set (rep1)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
#   add_row(Method = "STAR (encode)", Graph = "1000g (no-CEU)", Reads = "Training set (rep1)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0) %>%
#   add_row(Method = "STAR (encode)", Graph = "1000g (no-CEU)", Reads = "Training set (rep2)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ0_frac", Frac = 0) %>%
#   add_row(Method = "STAR (encode)", Graph = "1000g (no-CEU)", Reads = "Training set (rep2)", count = 0, MapQ0 = 0, MapQ1 = 0, Filter = "MapQ1_frac", Frac = 0)
#   
# overlap_data_mir_mapq_bar$Filter = recode_factor(overlap_data_mir_mapq_bar$Filter, 
#                                                    "MapQ0_frac" = "All",
#                                                    "MapQ1_frac" = "MapQ > 0")
# 
# overlap_data_mir_mapq_bar$Graph = recode_factor(overlap_data_mir_mapq_bar$Graph, 
#                                                   "Linear" = "Linear",
#                                                   "1000g (no-CEU)" = "1000g\n(no-CEU)")
# 
# pdf("plots/micro_rna/real_overlap_mir_mapq_bar.pdf", height = 4, width = 5, pointsize = 12)
# ggplot() +
#   geom_bar(data = overlap_data_mir_mapq_bar[overlap_data_mir_mapq_bar$Filter == "MapQ > 0",], aes(Graph, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge()) +
#   geom_bar(data = overlap_data_mir_mapq_bar[overlap_data_mir_mapq_bar$Filter == "All",], aes(Graph, y = Frac, fill = Method, alpha = Filter), stat = "identity", width = 0.5, position = position_dodge(), alpha = 0.5) +
#   scale_fill_manual(values = wes_cols_mir) +
#   scale_alpha_manual(name = "Filter", values = c(0.5, 1), labels = c("Unfiltered", "MapQ > 0"), drop = F) + 
#   scale_y_continuous(limits=c(0, 1), oob = rescale_none) +
#   #facet_grid(cols = vars(Reads)) +
#   xlab("") +
#   ylab("Mapping rate") +
#   theme_bw() +
#   theme(strip.background = element_blank()) +
#   theme(text = element_text(size=12))
# dev.off()
# 
# 
