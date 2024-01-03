
plot_data_column = function (data, column,size, colors ,sample ) {
  ggplot(data=data, aes_string(x="sites", y=column,color=sample)) +
    geom_point(size = size) +
    scale_color_manual(values=colors) +
    labs(title = column) +
    ylab(column) +
    theme_minimal() #+
    #ylim(range(data[["Y"]]))
}


plot_proc_ggplot = function(plot_Y, names_proc, size = 0.5, sample = "sample"){
  n = length(unique(plot_Y[[sample]]))
  colors = rep(brewer.pal(8,"Dark2"),n)[1:n]
  
  fig = lapply(names_proc,plot_data_column, data = plot_Y, size = size, colors, sample)
  ggarrange(plotlist=fig,ncol=1,common.legend = TRUE, legend="bottom")
}
