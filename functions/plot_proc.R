


plot_proc = function(plot_Y, names_proc,size=2){
  n = length(levels(plot_Y[["sample"]]))
  colors = rep(brewer.pal(8,"Dark2"),n)[1:n]
  fig = list()
  length(fig) = length(names_proc)
  fig[[1]] = plot_Y%>%
    group_by(colors)%>%
    plot_ly(x=~sites, color= ~colors, legendgroup=~colors)%>%
    add_trace(y= plot_Y[[names_proc[1]]], mode = "markers", type = "scatter", marker = list(size=size), name = ~sample)%>%#, marker = list(size=size)) %>%
    layout(
      annotations = list(
        text = names_proc[1],
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
      ))
  
  if(length(names_proc)>1){
    for ( i in 2:length(names_proc)){
      fig[[i]] = plot_Y%>%
        group_by(colors)%>%
        plot_ly(x=~sites, color= ~colors, legendgroup=~colors)%>%
        add_trace(y= plot_Y[[names_proc[i]]], mode = "markers", type = "scatter", showlegend=F, marker = list(size=size), name = ~sample) %>%
        layout(
          annotations = list(
            text = names_proc[i],
            xref = "paper",
            yref = "paper",
            yanchor = "bottom",
            xanchor = "center",
            align = "center",
            x = 0.5,
            y = 1,
            showarrow = FALSE
          ))
    }
  }
  
 
  
  subplot(fig,nrows = length(fig), shareX = T,titleY = TRUE)
}
