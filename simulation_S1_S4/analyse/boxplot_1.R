boxplot_1 = function(value,
                     ylab = "RMSE",

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
                     title = "Root Mean Square Error",

                     col_boxplot = "black",
                     fill_box = "grey",
                     transparence_box = 0.8,

                     y_min = 0,
                     y_max = 1
                     ) {
  if (ncol(value) > 2) {
    stop("More than 2 columns in value")
  } else{
    df = data.table(value)
    n = nrow(df)
    lab = c(rep(colnames(df)[1], n), rep(colnames(df)[2], n))
    a = df[[1]]
    b = df[[2]]
    df = data.table(col = c(a, b), label = factor(lab,levels = unique(lab)))


    p = ggplot(df, aes(x = label, y = col, fill = label)) +
      geom_boxplot(color=col_boxplot, alpha=transparence_box) + xlab("") + ylab(ylab) + ggtitle(title) +
      theme_minimal() +
      scale_fill_manual(values=fill_box) +
      ylim(y_min, y_max)+
      theme(
        text = element_text(size=textSize),
        legend.position="none",
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
      ggsave(path, width = width, height = height,device=cairo_ps )
    } else{
      return(p)
    }

  }
}
