barplot_1 = function(tab_time,
                     
                     save_pdf = F,
                     path = "",
                     width = 6,
                     height = 5,
                     
                     colorTitle = "black",
                     sizeTitle = 15,
                     formeTitle = "bold.italic",
                     colorAxe = "black",
                     sizeAxe = 10,
                     formeAxe = "bold",
                     textSize = 15,
                     ylab = "Time in sec",
                     title = "Total of computation times",
                     
                     col_boxplot = "black",
                     fill_box = "grey",
                     transparence_box = 0.8
) {
    
    p = ggplot() + geom_col(data = tab_time,
                            aes(x = simulations, y = Time, fill = Method),
                            position = "dodge") +
      scale_fill_manual(values=fill_box) +
      labs(title = title,
           x = "Simulation",
           y = ylab) +
      theme_minimal()  +
      theme(
        text = element_text(size=textSize),
        plot.title = element_text(
          hjust = 0.5,
          color = colorTitle,
          size = sizeTitle,
          face = formeTitle
        ),
        axis.title.x = element_text(
          color = colorAxe,
          size = sizeAxe,
          face = formeAxe
        ),
        axis.title.y = element_text(
          color = colorAxe,
          size = sizeAxe,
          face = "bold"
        )
      )
    if (save_pdf) {
      p
      ggsave(path, width = width, height = height, )
    } else{
      return(p)
    }
    
}