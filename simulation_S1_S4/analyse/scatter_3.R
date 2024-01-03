scatter_3 = function(sites,
                     Y,
                     Y_obs,
                     proc_a,
                     proc_b,
                     x,

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
                     Title = "Two samples",

                     low = "#349be8",
                     high = "#cc0000",
                     point_size = 2){


  df = data.table(sites = sites,
                  Y_A = Y[proc_a, ],
                  Y_A_obs = Y_obs[proc_a, ],
                  Y_B = Y[proc_b, ],
                  Y_B_obs = Y_obs[proc_b, ])


  df = data.table(df %>% pivot_longer(!sites, names_to = "origine", values_to = "values"))

  df = data.table(df, Type = factor(rep(c("Obs","True"),nrow(df)/2),labels = c("True","Obs")),
                  X = rep(rep(c(x[proc_a],x[proc_b]),each=2),nrow(df)/4))

 p = ggplot(df, aes(
    x = sites,
    y = values,
    color = X,
    shape = Type,
    size = Type
  )) +
   geom_point(alpha = 1) +
   scale_color_gradient(low = low, high = high, limits = c(0,1))+
   scale_shape_manual(values = c(1, 15)) +
   scale_size_manual(guide="none", values = c(1,5)) +
   guides(shape = guide_legend(override.aes = list(size=point_size))) +
    labs(title = Title,
         x = "Sites s",
         y = "Y(s)")+
    theme_minimal() +
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
