
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
    group_by(Method, Graph, Threshold, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    arrange(desc(MapQ), .by_group = T) %>%
    mutate(TPcs = cumsum(TP), FPcs = cumsum(FP)) %>%
    mutate(N = max(TPcs) + max(FPcs)) %>% 
    mutate(Sensitivity = (FPcs + TPcs) / N, Precision = TPcs / (FPcs + TPcs)) 
  
  overlap_data_roc %>%
    filter(MapQ == 30) %>%
    print()
  
  min_lim_xy <- min(c(min(overlap_data_roc$Sensitivity), min(overlap_data_roc$Precision)))
  
  p <- overlap_data_roc %>%
    ggplot(aes(y = Precision, x = Sensitivity, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    facet_grid(cols = vars(Threshold)) +
    scale_color_manual(values = cols) +
    coord_fixed() +
    xlim(c(min_lim_xy, 1)) +
    ylim(c(min_lim_xy, 1)) +
    xlab("Mapping sensitivity") +
    ylab("Mapping precision") +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=18))
  print(p)
  
  min_lim_x <- min(overlap_data_roc$Sensitivity)

  p <- overlap_data_roc %>%
    ggplot(aes(y = log10(1 - Precision), x = Sensitivity, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    facet_grid(cols = vars(Threshold)) +
    scale_color_manual(values = cols) +
    xlim(c(min_lim_x, 1)) +
    xlab("Mapping sensitivity") +
    ylab("Mapping error (log10)") +
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=18)) 
  print(p)
}

plotMeanOverlap <- function(overlap_data, cols) {
  
  overlap_data_mean <- overlap_data %>% 
    mutate(OverlapCount = Overlap * Count) %>%
    group_by(Method, Graph, Threshold, MapQ) %>%
    summarise(SumCount = sum(Count), SumOverlapCount = sum(OverlapCount)) %>% 
    mutate(MeanOverlap = SumOverlapCount/ SumCount)
  
  p <- overlap_data_mean %>%
    ggplot(aes(y = MeanOverlap, x = MapQ, color = Method, linetype = Graph, shape = Graph)) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    scale_color_manual(values = cols) +
    ylab("Mean overlap") +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=18))
  print(p)
}

plotMapQ <- function(overlap_data, cols) {

  overlap_data_mapq <- overlap_data %>% 
    group_by(Method, Graph, Threshold, MapQ) %>%
    summarise(TP = sum(TP), FP = sum(FP)) %>% 
    mutate(LogMapQError = -1 * MapQ / 10) %>%
    mutate(LogEstError = log10(1 - TP / (FP + TP)))
  
  p <- overlap_data_mapq %>%
    ggplot(aes(y = LogEstError, x = LogMapQError, color = Method, linetype = Graph, shape = Graph)) +
    geom_abline(intercept = 0) +
    geom_line(size = 1) + 
    geom_point(size = 1.5) +
    facet_grid(cols = vars(Threshold)) +
    scale_color_manual(values = cols) +
    ylab("Log estimated error") +
    xlab("Log mapQ error") +
    xlim(c(-6, 0)) +
    ylim(c(-6, 0)) +
    coord_fixed() +
    theme_bw() +
    theme(strip.background = element_blank()) +
    theme(text = element_text(size=18))
  print(p)
}

plotOverlapBenchmark <- function(overlap_data, cols, filename) {

  pdf(filename, width = 9)
  
  overlap_data <- overlap_data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
  plotRocCurve(overlap_data, cols)
  plotMapQ(overlap_data, cols)
  
  dev.off() 
}

plotDistanceBenchmark <- function(distance_data, cols, filename) {
  
  pdf(filename, width = 9)
  
  distance_data <- distance_data %>%
    mutate(TP = Count * Correct) %>% 
    mutate(FP = Count * !Correct)
  
  plotRocCurve(distance_data, cols)
  plotMapQ(distance_data, cols)
  
  dev.off() 
}
