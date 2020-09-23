
# utils.R

library("ggrepel")
library("scales")

wes_cols <- c(rev(wes_palette("Rushmore1")[c(2,1,3,4,5)]), wes_palette("Zissou1")[c(1)])

num_reads <- list()
num_reads["470"] <- 115359773
num_reads["SRR1153470"] <- 115359773

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

plotRocCurve <- function(data, cols, log = T) {
  
  set.seed(1234)
  
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
  
  if (log) {
    
    a <- annotation_logticks(sides = "l")
    a$data <- data.frame(x = NA, FacetCol = c(data_roc$FacetCol[1]))
    
    p <- data_roc %>%
      ggplot(aes(y = -1 * log10(1 - Precision), x = Sensitivity, color = Method, linetype = Graph, shape = Graph, label = MapQ)) +
      a +
      geom_line(size = 0.75) +
      geom_point(data = subset(data_roc, MapQ == 0 | MapQ == 1 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60 | MapQ == 255), size = 1.75, alpha = 1) +
      geom_text_repel(data = subset(data_roc, MapQ == 0 | MapQ == 1 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60| MapQ == 255), size = 3, fontface = 2) +
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
  
  } else {
  
    p <- data_roc %>%
      ggplot(aes(y = Precision, x = Sensitivity, color = Method, linetype = Graph, shape = Graph, label = MapQ)) +
      geom_line(size = 0.75) +
      geom_point(data = subset(data_roc, MapQ == 0 | MapQ == 1 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60 | MapQ == 255), size = 1.75, alpha = 1) +
      geom_text_repel(data = subset(data_roc, MapQ == 0 | MapQ == 1 | (MapQ == 42 & grepl("Bowtie2", Method)) | MapQ == 60| MapQ == 255), size = 3, fontface = 2) +
      facet_grid(FacetRow ~ FacetCol) +
      scale_color_manual(values = cols) +
      xlim(c(min_lim_x, 1)) +
      xlab("Mapping sensitivity") +
      ylab("Fraction mapped reads overlapping microRNAs") +
      theme_bw() +
      theme(aspect.ratio=1) +
      theme(strip.background = element_blank()) +
      theme(text = element_text(size=14))
    print(p)
  }
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

plotOverlapBenchmark <- function(data, cols, filename, log = T) {

  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(data, cols)
#  dev.off() 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  plotRocCurve(data, cols, log)
  dev.off() 
}

plotDistanceBenchmark <- function(data, cols, filename) {
  
  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(data, cols)
#  dev.off() 
  
  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  plotRocCurve(data, cols)
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
