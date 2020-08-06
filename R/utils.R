
# utils.R

printHeader <- function() {

  args <- commandArgs()
  cript_dir <- dirname(sub("--file=", "", args[4]))
  print(script_dir)
  
  print(args)
  system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
  system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))
}

plotRocCurve <- function(data, cols) {
  
  data_roc <- data %>% 
    group_by(Reads, Method, Graph, Threshold, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
    mutate(N = max(TPcs) + max(FPcs)) %>% 
    mutate(Sensitivity = (FPcs + TPcs) / N, Precision = TPcs / (FPcs + TPcs)) 
  
  data_roc %>%
    filter(MapQ >= 60) %>%
    print()
  
  min_lim_xy <- min(c(min(data_roc$Sensitivity), min(data_roc$Precision)))
  
  # p <- data_roc %>%
  #   ggplot(aes(y = Precision, x = Sensitivity, color = Method, linetype = Graph, shape = Graph)) +
  #   geom_line(size = 1) + 
  #   geom_point(size = 1.5) +
  #   facet_grid(Reads ~ Threshold) +
  #   scale_color_manual(values = cols) +
  #   coord_fixed() +
  #   xlim(c(min_lim_xy, 1)) +
  #   ylim(c(min_lim_xy, 1)) +
  #   xlab("Mapping sensitivity") +
  #   ylab("Mapping precision") +
  #   theme_bw() +
  #   theme(strip.background = element_blank()) +
  #   theme(text = element_text(size=14))
  # print(p)
  
  min_lim_x <- min(data_roc$Sensitivity)
  min_lim_y <- min(log10(1 - data_roc$Precision))
  
  p <- data_roc %>%
    ggplot(aes(y = log10(1 - Precision), x = Sensitivity, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    facet_grid(Threshold ~ Reads) +
    scale_color_manual(values = cols) +
    xlim(c(min_lim_x, 1)) +
    #ylim(c(min_lim_y, 0)) +
    xlab("Mapping sensitivity") +
    ylab("Mapping error (log10)") +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14)) 
  print(p)
}

plotMapQ <- function(data, cols) {

  data_mapq <- data %>% 
    group_by(Reads, Method, Graph, MapQ) %>%
    summarise(Count = sum(Count)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(Countcs = cumsum(Count)) %>%
    mutate(N = max(Countcs)) %>% 
    mutate(Sensitivity = Countcs / N) 
  
  p <- data_mapq %>%
    ggplot(aes(y = Sensitivity, x = MapQ, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    facet_grid(cols = vars(Reads)) +
    scale_color_manual(values = cols) +
    xlab("Mapping quality") +
    ylab("Mapping sensitivity") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14))
  print(p)
}

plotOverlapBenchmark <- function(data, cols, filename) {

  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(data, cols)
#  dev.off() 
  
  png(paste(filename, ".png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  plotRocCurve(data, cols)
  dev.off() 
}

plotDistanceBenchmark <- function(data, cols, filename) {
  
  data <- data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(data, cols)
#  dev.off() 
  
  png(paste(filename, ".png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  plotRocCurve(data, cols)
  dev.off() 
}

plotMapQBenchmark <- function(data, cols, filename) {
  
  #  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
  #  plotRocCurve(data, cols)
  #  dev.off() 
  
  png(paste(filename, ".png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  plotMapQ(data, cols)
  dev.off() 
}
