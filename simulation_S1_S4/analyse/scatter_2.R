scatter_2 = function(sites,
                     ind_na,
                     Y,
                     Y_GaSP,
                     Y_OLS_GaSP,
                     proc,

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
                     Title = "Residuals",
                     
                     col_gasp = "#f0d8c9",
                     col_ols_gasp = "#f46d61",
                     point_size = 2,
                     point_shape = 16){
  
  
  df = data.table(sites = sites[ind_na],
                  GaSP = (Y-Y_GaSP)[proc, ind_na],
                  OLS_GaSP = (Y-Y_OLS_GaSP)[proc, ind_na])
  
  df = data.table(df %>% pivot_longer(!sites, names_to = "origine", values_to = "values"))
  df$origine = factor(df$origine, levels = c("OLS_GaSP","GaSP"))
  
  
  
  
  p = ggplot(df, aes(x=sites, y=values, color = origine))+ 
    geom_point(shape=point_shape,
               alpha=1,
               size=point_size
               )+
    scale_color_manual(paste0("sample_",proc), values=c(col_gasp, col_ols_gasp))+
    geom_hline(yintercept = 0, color = "black") +
    labs(title = Title,
         x = "Sites",
         y = "Residuals")+
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
