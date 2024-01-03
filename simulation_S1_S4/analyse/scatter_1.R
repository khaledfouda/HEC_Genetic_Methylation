scatter_1 = function(sites,
                     Y,
                     Y_obs,
                     proc_a,
                     proc_b,
                     x = x,
                     
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
                     
                     col_y_a = "#2986cc",
                     col_y_a_obs = "#d0e6f6",
                     col_y_b = "#eea92f",
                     col_y_b_obs = "#f6e8d0",
                     point_size = 2){

    
    df = data.table(sites = sites,
                    Y_A = Y[proc_a, ],
                    Y_A_obs = Y_obs[proc_a, ],
                    Y_B = Y[proc_b, ],
                    Y_B_obs = Y_obs[proc_b, ])
    
    
    df = data.table(df %>% pivot_longer(!sites, names_to = "origine", values_to = "values"))
    

    p = ggplot() +
      geom_point(
        data = subset(df, grepl("^Y_A", origine)),
        aes(
          x = sites,
          y = values,
          shape = origine,
          color = origine
        ),
        size = point_size
      ) +
      scale_shape_manual(
        name = paste0("sample_", proc_a),
        values = c(Y_A = 1, Y_A_obs = 16),
        labels = c(Y_A = "y_true", Y_A_obs = "y_pred"),
        guide = guide_legend(order = 1)
      ) +
      scale_color_manual(
        name = paste0("sample_", proc_a),
        values = c(Y_A = col_y_a, Y_A_obs = col_y_a_obs),
        labels = c(Y_A = "y_true", Y_A_obs = "y_pred"),
        guide = guide_legend(order = 1)
      ) +
      new_scale_color() +
      new_scale("shape") +
      geom_point(
        data = subset(df, grepl("^Y_B", origine)),
        aes(
          x = sites,
          y = values,
          shape = origine,
          color = origine
        ),
        size = point_size
      ) +
      scale_shape_manual(
        name = paste0("sample_", proc_b),
        values = c(Y_B = 1, Y_B_obs = 16),
        labels = c(Y_B = "y_true", Y_B_obs = "y_pred"),
        guide = guide_legend(order = 2)
      ) +
      scale_color_manual(
        name = paste0("sample_", proc_b),
        values = c(Y_B = col_y_b, Y_B_obs = col_y_b_obs),
        labels = c(Y_B = "y_true", Y_B_obs = "y_pred"),
        guide = guide_legend(order = 2)
      )+
      labs(title = Title,
           x = "Sites",
           y = "Gaussian processes")+
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