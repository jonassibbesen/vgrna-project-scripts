
# utils.R

library("ggrepel")
library("scales")

wes_cols <- c(rev(wes_palette("Rushmore1")[c(2,1,3,4,5)]), wes_palette("Zissou1")[c(1)])

data_set1 <- "SRR1153470_uni"
data_set2 <- "470"
data_set3 <- "SRR1153470"

printHeader <- function() {

  args <- commandArgs()
  cript_dir <- dirname(sub("--file=", "", args[4]))
  print(script_dir)
  
  print(args)
  system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
  system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))
}

plotRocCurveMapq <- function(data, cols) {
  
  set.seed(1234)
  
  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
  data_roc <- data %>% 
    mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
    group_by(Method, Graph, FacetRow, FacetCol, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
    mutate(N = max(TPcs) + max(FPcs)) %>% 
    mutate(Sensitivity = (FPcs + TPcs) / N, Precision = TPcs / (FPcs + TPcs)) %>%
    filter(MapQ >= 0)
  
  min_lim_x <- min(data_roc$Sensitivity)
  
  a <- annotation_logticks(sides = "l")
  a$data <- data.frame(x = NA, FacetCol = c(as.character(data_roc$FacetCol[1])))
  
  data_roc %>% filter(MapQ >= 60) %>% print(n = 100)
  
  p <- data_roc %>%
    ggplot(aes(y = -1 * log10(1 - Precision), x = Sensitivity, color = Method, linetype = Graph, shape = Graph, label = MapQ)) +
    a +
    geom_line(size = 1) +
    geom_point(data = subset(data_roc, MapQ == 0 | MapQ == 1 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60 | MapQ == 255), size = 2, alpha = 1) +
    geom_text_repel(data = subset(data_roc, MapQ == 0 | MapQ == 1 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60| MapQ == 255), size = 3.5, fontface = 2) +
    scale_y_continuous(breaks = seq(1, 4), labels = c(0.9, 0.99, 0.999, 0.9999)) + 
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    xlim(c(min_lim_x, 1)) +
    xlab("Mapping sensitivity") +
    ylab("Mapping accuracy") +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14)) 
  print(p)   
}

plotRocCurveOvl <- function(data, cols) {
  
  set.seed(1234)
  
  data_roc <- list()
  
  for (i in seq(0.05, 1, 0.05)) {
    
    data_roc_thres <- data %>% 
      filter(IsMapped > 0) %>% 
      add_column(Threshold = i) %>%
      mutate(Correct = Overlap > i) %>%
      mutate(TP = Count * Correct) %>% 
      mutate(FP = Count * !Correct) %>% 
      group_by(Method, Graph, FacetRow, FacetCol, Threshold) %>%
      summarise(TP = sum(TP), FP = sum(FP)) %>% 
      arrange(Threshold, .by_group = T) %>%
      mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
      mutate(Precision = TPcs / (FPcs + TPcs)) 
    
    data_roc[[as.character(i)]] <- data_roc_thres
  }
  
  data_roc <- do.call(rbind, data_roc)
  
  a <- annotation_logticks(sides = "l")
  a$data <- data.frame(x = NA, FacetCol = c(as.character(data_roc$FacetCol[1])))
  
  print(head(data_roc))
  print(c(as.character(data_roc$FacetCol[1])))
  print(a$data)
  
  p <- data_roc %>%
    ggplot(aes(y = -1 * log10(1 - Precision), x = Threshold, color = Method, linetype = Graph, shape = Graph, label = Threshold)) +
    a +
    geom_line(size = 0.75) +
    geom_text_repel(data = subset(data_roc, Threshold == 0.5), size = 3, fontface = 2) +
    scale_y_continuous(breaks = seq(1, 4), labels = c(0.9, 0.99, 0.999, 0.9999)) + 
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    xlim(c(0, 1)) +
    xlab("Overlap threshold") +
    ylab("Mapping accuracy") +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14)) 
  print(p)   
}

plotF1Curve <- function(data, cols, log = T) {
  
  data_f1 <- data %>% 
    mutate(MapQ = ifelse(MapQ > 60, 60, MapQ)) %>% 
    mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
    group_by(Method, Graph, FacetRow, FacetCol, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
    mutate(N = max(TPcs) + max(FPcs)) %>% 
    mutate(Sens = (FPcs + TPcs) / N, Pre = TPcs / (FPcs + TPcs)) %>%
    mutate(F1 = (2 * Sens * Pre) / (Sens + Pre)) %>%
    filter(MapQ >= 0)
  
  p <- data_f1 %>%
    ggplot(aes(y = F1, x = MapQ, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 0.75, alpha = 0.90) +
    geom_point(size = 1.25, alpha = 0.90) +
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    xlab("Mapping quality") +
    ylab("Mapping F1 score") +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=18))
  print(p)
}

plotMapQCurve <- function(data, cols) {

  data_mapq <- data %>% 
    mutate(MapQ = ifelse(MapQ > 60, 60, MapQ)) %>% 
    mutate(MapQ = ifelse(IsMapped, MapQ, -1)) %>% 
    group_by(Method, Graph, FacetRow, FacetCol, MapQ) %>%
    summarise(Count = sum(Count)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(Countcs = cumsum(Count)) %>%
    mutate(N = max(Countcs)) %>% 
    mutate(Sensitivity = Countcs / N) %>%
    filter(MapQ >= 0) 
  
  p <- data_mapq %>%
    ggplot(aes(y = Sensitivity, x = MapQ, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 0.75, alpha = 0.90) +
    geom_point(size = 1.25, alpha = 0.90) +
    facet_grid(FacetRow ~ FacetCol) +
    scale_color_manual(values = cols) +
    xlab("Mapping quality") +
    ylab("Mapping sensitivity") +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=18))
  print(p)
}

plotOverlapBenchmarkMapQ <- function(data, cols, filename) {

#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(data, cols)
#  dev.off() 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  plotRocCurveMapq(data, cols)
  dev.off() 
}

plotOverlapBenchmarkOvl <- function(data, cols, filename) {
  
  #  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  #  plotRocCurve(data, cols)
  #  dev.off() 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  plotRocCurveOvl(data, cols)
  dev.off() 
}

plotDistanceBenchmarkMapQ <- function(data, cols, filename) {
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(data, cols)
#  dev.off() 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  plotRocCurveMapq(data, cols)
  dev.off() 
}

plotMapQBenchmark <- function(data, cols, filename) {
  
  #  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  #  plotRocCurve(data, cols)
  #  dev.off() 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 7, pointsize = 12)
  plotMapQCurve(data, cols)
  dev.off() 
}
