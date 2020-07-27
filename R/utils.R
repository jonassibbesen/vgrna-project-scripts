
# utils.R

printHeader <- function() {

  args <- commandArgs()
  cript_dir <- dirname(sub("--file=", "", args[4]))
  print(script_dir)
  
  print(args)
  system(paste(c("git", "-C", script_dir, "rev-parse", "HEAD"), collapse = " "))
  system(paste(c("git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"), collapse = " "))
}

plotRocCurve <- function(overlap_data, cols) {
  
  overlap_data_roc <- overlap_data %>% 
    group_by(Reads, Method, Graph, Threshold, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
    mutate(N = max(TPcs) + max(FPcs)) %>% 
    mutate(Sensitivity = (FPcs + TPcs) / N, Precision = TPcs / (FPcs + TPcs)) 
  
  overlap_data_roc %>%
    filter(MapQ >= 60) %>%
    print()
  
  min_lim_xy <- min(c(min(overlap_data_roc$Sensitivity), min(overlap_data_roc$Precision)))
  
  # p <- overlap_data_roc %>%
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
  
  min_lim_x <- min(overlap_data_roc$Sensitivity)
  min_lim_y <- min(log10(1 - overlap_data_roc$Precision))
  
  p <- overlap_data_roc %>%
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

plotMapQ <- function(overlap_data, cols) {

  overlap_data_mapq <- overlap_data %>% 
    group_by(Reads, Method, Graph, Threshold, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    mutate(LogMapQError = -1 * MapQ / 10) %>%
    mutate(LogEstError = log10(1 - TP / (FP + TP)))
  
  p <- overlap_data_mapq %>%
    ggplot(aes(y = LogEstError, x = LogMapQError, color = Method, linetype = Graph, shape = Graph)) +
    geom_abline(intercept = 0) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    facet_grid(Threshold ~ Reads) +
    scale_color_manual(values = cols) +
    ylab("Log estimated error") +
    xlab("Log mapQ error") +
    xlim(c(-6, 0)) +
    ylim(c(-6, 0)) +
    coord_fixed() +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=14))
  print(p)
}

plotOverlapBenchmark <- function(overlap_data, cols, filename) {

  overlap_data <- overlap_data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(overlap_data, cols)
#  dev.off() 
  
  png(paste(filename, ".png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  plotRocCurve(overlap_data, cols)
  dev.off() 
}

plotDistanceBenchmark <- function(distance_data, cols, filename) {
  
  distance_data <- distance_data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
#  pdf(paste(filename, ".pdf", sep = ""), height = 6, width = 9, pointsize = 12)
#  plotRocCurve(distance_data, cols)
#  dev.off() 
  
  png(paste(filename, ".png", sep = ""), height = 6, width = 9, units = "in", pointsize = 12, res = 300)
  plotRocCurve(distance_data, cols)
  dev.off() 
}
