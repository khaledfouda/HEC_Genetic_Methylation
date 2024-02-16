compute_p_values_and_plot_multivar <- function(chromosome, load_p_values=TRUE, Male.Only=TRUE, plot = TRUE,
                              alpha=0.05, min_freq=1, middle_point=FALSE, floor_by=1e7, note = "",
                              correction = function(alpha, N) (alpha / N)){

   #' This function either computes the p_values as discussed in the document or plot them or both.
   #' 
   Y = readRDS(paste0("new_data/Ydat_common_",chromosome, note, ".rds"))
   X = readRDS(paste0("new_data/Xdat_common_",chromosome, note, ".rds")) %>% as.data.frame() #%>% 
   sites = readRDS(paste0("new_data/sites_common_",chromosome, note, ".rds"))
   
   
   if(Male.Only == TRUE){
      male_indices = which(X$MALE == 1)
      Y = Y[male_indices,]
      X = X[male_indices,]
   }
   
   site_order = order(sites) 
   Y = Y[,site_order]
   sites = sites[site_order]
   
   dim(Y)
   N = ncol(Y)
   K = nrow(Y)
   #--------------------------------
   if(load_p_values == TRUE){
      p_values = readRDS(paste0("new_data/p_values_Mutlivar_",chromosome, note, ".rds"))
   }else{
      library(parallel)
      no_cores <- 12             
      cl <- makeCluster(no_cores)
      
      Age = X$AGE
      Sex = X$MALE
      #clusterExport(cl, varlist = c("Y"))
      
      p_values <- parSapply(cl, 1:ncol(Y), function(i) {
         model <- lm(Y[, i] ~ Age + Sex) 
         c(Age = summary(model)$coefficients["Age", "Pr(>|t|)"], 
           Sex = summary(model)$coefficients["Sex", "Pr(>|t|)"])
      }) 
      
      p_values <- as.data.frame(t(p_values), stringsAsFactors = FALSE)
      colnames(p_values) <- c("Age_p_value", "Sex_p_value")
      
      stopCluster(cl)
      saveRDS(p_values, paste0("new_data/p_values_Mutlivar_",chromosome, note, ".rds"))
   }
   
   if(plot == TRUE){
      vars = c("Age", "Sex")
      for(i in 1:2){
         filename = paste0("./graphs/manhattan_plot_Multivar_",i,"_", chromosome,"_alpha_",alpha, note, ".png")
         print(paste("Saving plot to ",filename))
         corrected_alpha = correction(alpha, N)
         
         dmr_regions = get_dmr_regions(p_values[,i], sites, alpha= corrected_alpha, min_freq = min_freq, return_seq = T,
                                       middle_point = middle_point,floor_by = floor_by)
         
         
         print(length(dmr_regions))
         data.frame(x = 1:N, Site= sites, NegLogP = - log10(p_values[,i])) %>% 
            mutate(color = ifelse(Site %in% dmr_regions,"Correlated","Uncorrelated")) ->
            manh.dat
         print(sum(manh.dat$color == "Correlated"))
         #---------------------------------------------------------------------
         #--------------------------------------------------------------------
         gtitle <- paste0("Chromosome ", strsplit(chromosome,"chr")[[1]][2], "; alpha = ",
                          alpha, "; corrected alpha = ", corrected_alpha, "; ",vars[i])
         manh.dat %>% 
            ggplot(aes(Site, NegLogP, color=color, fill=color)) +  
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c("Correlated" = "blue", "Uncorrelated" = "black")) +
            theme_minimal() +
            #scale_x_continuous(limits=c(0,2e8))+
            xlab("Methylation Site") +
            ylab("-log10(P-value)") +
            geom_hline(yintercept = -log10(corrected_alpha), linetype = "dashed", color = "red")+
            geom_hline(yintercept = -log10((.05/ N)), linetype = "dashed", color = "red")+ 
            geom_hline(yintercept = -log10((0.1/ N)), linetype = "dashed", color = "red")+ 
            geom_hline(yintercept = -log10((0.2/ N)), linetype = "dashed", color = "red")+ 
            ggtitle(gtitle, subtitle=note) -> p 
         
         ggsave(filename = filename,
                plot = p, width = 10, height = 6, dpi = 300) 
         print("Saved.")
      }
   }
}

#------------------------
# Example use:
# for(chr in paste0("chr",c(7,8,11,12,17))){
#    for(alpha in c(0.05, .1, .2))
#       compute_p_values_and_plot(chr, load_p_values = TRUE, alpha=alpha, min_freq=1, middle_point=TRUE, floor_by=1e6, plot=TRUE)
#    
# }
#--------------------------------------
