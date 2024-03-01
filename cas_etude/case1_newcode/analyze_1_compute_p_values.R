compute_p_values_case1 <-
 function(Y,
          sites,
          X,
          load_p_values = TRUE,
          plot = TRUE,
          alpha = 1e-4,
          note = "",
          normalize = TRUE,
          floor_by = 1e3,
          correction=FALSE) {
  dim(Y)
  N = ncol(Y)
  K = nrow(Y)
  filename = paste0("new_data/case1_p_values", note, ".rds")
  if (normalize)
   X = (X - mean(X)) / sd(X)
  
  if(correction) alpha = alpha / N
  
  if (load_p_values) {
   p_values = readRDS(filename)
  } else{
   require(parallel)
   no_cores = 12
   cl = makeCluster(no_cores)
   
   p_values <- parSapply(cl, 1:ncol(Y), function(i) {
    summary(lm(Y[, i] ~ X))$coefficients["X", "Pr(>|t|)"]
   })
   stopCluster(cl)
   saveRDS(p_values, filename)
   print(paste0("p_values saved to ", filename))
  }
  if (!plot)
   return(p_values)
  
  plot_file = paste0("graphs/case1_Manhattan_plot", note, ".png")
  print(paste0("Saving plot to ", plot_file))
  
  dmr_regions <- get_dmr_regions_case1(p_values,
                                       sites,
                                       alpha,
                                       floor_by)
  
  print(paste0("Length of DMR = ", length(dmr_regions)))
  manh.dat <- data.frame(x = 1:N,
                         site = sites,
                         NegLogP = -log10(p_values)) %>%
   mutate(color = ifelse(site %in% dmr_regions, "Significant", "Not Significant"))
  
  print(paste0(
   "Number of significant sites is ",
   sum(manh.dat$color == "Significant")
  ))
  
  gtitle <- paste0("Alpha = ", alpha, " - case 1")
  
  manh.dat %>%
   ggplot(aes(site, NegLogP, color = color, fill = color)) +
   geom_point(alpha = 0.6) +
   scale_color_manual(values = c(
    "Significant" = "blue",
    "Not Significant" = "black"
   )) +
   theme_minimal() +
   xlab("Methylation Site") +
   ylab("-log10(p-value)") +
   geom_hline(
    yintercept = -log10((.05 / N)),
    linetype = "dashed",
    color = "red"
   ) +
   geom_hline(
    yintercept = -log10((0.1 / N)),
    linetype = "dashed",
    color = "red"
   ) +
   geom_hline(
    yintercept = -log10((0.2 / N)),
    linetype = "dashed",
    color = "red"
   ) +
   geom_hline(
    yintercept = -log10(alpha),
    linetype = "dashed",
    color = "darkgreen"
   ) +
   ggtitle(gtitle, subtitle = note) -> p
  
  ggsave(
   filename = plot_file,
   plot = p,
   width = 10,
   height = 6,
   dpi = 300
  )
  
  print("Plot saved.")
  return(manh.dat)
 }
