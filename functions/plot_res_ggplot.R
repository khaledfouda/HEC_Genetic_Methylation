plot_data_pred = function (data, column, colors, sample) {
  ggplot(data=data, aes_string(x="Y", y=column,color=sample)) +
    geom_point() +
    scale_color_manual(values=colors) +
    labs(title = paste(column,": RMSE =" ,round(RMSE(data$Y,data[[column]]),3))) +
    ylab(column) +
    geom_abline(intercept = 0, slope = 1, col = "red") + 
    theme_minimal() +
    theme(plot.title = element_text(size=10))
  
}

plot_res_ggplot = function(plot_Y, proc,names_proc = proc, sample = "sample"){
  n = length(unique(plot_Y[[sample]]))
  colors = rep(brewer.pal(12,"Paired"),n)[1:n]
  fig = lapply(proc,plot_data_pred, data = plot_Y[is.na(Y_obs)],colors = colors, sample = sample)
  ggarrange(plotlist=fig, ncol=2,nrow=ceiling(length(fig)/2),common.legend = TRUE, legend="bottom")
  
}