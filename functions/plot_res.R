plot_res = function(plot_Y, proc,n_star, names_proc = proc){
  
  plot_Y_na = plot_Y[is.na(Y_obs)]
  fig = list()
  length(fig) = length(names_proc)
  
  titre = paste(names_proc[1],": R2 =" ,round(R2(plot_Y_na$Y,plot_Y_na[[proc[1]]]),3))
  fig[[1]] = plot_Y_na %>%
    group_by(sample) %>%
    plot_ly(x =  ~ Y,
            color = ~ sample,
            legendgroup =  ~ sample) %>%
    add_trace(y = plot_Y_na[[proc[1]]], mode = "markers",  showlegend=F,type = "scatter") %>%
    layout(
      annotations = list(
        text = titre,
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
      ),
      legend = list(title = list(text = "Sample")),
      xaxis = list(title = 'Y_true',
                   range = range(plot_Y_na$Y)),
      yaxis = list(title = 'Y_pred',
                   range = range(plot_Y_na$Y))) %>%
        add_trace(x = ~ Y, y = ~ Y, mode = "lines", type = "scatter", showlegend=F,line=list(color="red")) 
    
  
  for ( i in 2:length(names_proc)){
    titre = paste(names_proc[i],": R2 =" ,round(R2(plot_Y_na$Y,plot_Y_na[[proc[i]]]),3))
    fig[[i]] = plot_Y_na %>%
      group_by(sample) %>%
      plot_ly(x =  ~ Y,
              color = ~ sample,
              legendgroup =  ~ sample) %>%
      add_trace(y = plot_Y_na[[proc[i]]], mode = "markers", type = "scatter", showlegend=F) %>%
      layout(
        annotations = list(
          text = titre,
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        ),
        legend = list(title = list(text = "Sample")),
        xaxis = list(title = 'Y_true',
                     range = range(plot_Y_na$Y)),
        yaxis = list(title = 'Y_pred',
                     range = range(plot_Y_na$Y))) %>%
      add_trace(x = ~ Y, y = ~ Y, mode = "lines", type = "scatter", showlegend=F,line=list(color="red"))
  }
  
  subplot(fig,nrows = round(length(fig)/2),shareX = T, shareY = T)

}