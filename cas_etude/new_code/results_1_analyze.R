
get_results_data <- function(min_freq = 1)
{
   all_res <-data.frame()
   for(Age.Only in c(T,F)){
      for(alpha in c(.05,.1,.2)){
         print(paste("*_*_*_*_*_*_*_*_*_*_*_ Alpha = ",alpha,"*_*_*_*_*_*_*_*_*_*_*_"))
         for(chromosome in paste0("chr",c(7,8,11,12,17))){
            
            if(chromosome == "chr17" & Age.Only == FALSE & alpha == 0.1) next
            print("***********************************************************************")
            load(paste0("results/res_dmr_",chromosome,"_MaleOnly_TRUE_AgeOnly_",Age.Only ,"_",
                        alpha,"x",min_freq,".Rdata"))
            print(paste0("Chromosome ", strsplit(chromosome,"chr")[[1]][2], "; MALE ONLY;",
                         ifelse(Age.Only, " AGE ONLY", " ALL COVARIATES"), "; alpha = ",alpha,"; min_freq = ",min_freq," >>>>>>>>"))
            res %>%
               mutate(RMSE = round(RMSE, 3), R2 = round(R2, 3), RMSE_dmr = round(RMSE_dmr,3),
                      R2_dmr = round(R2_dmr, 3), time_min = round(time/60,1), time=NULL,
                      N = paste0(round(N/1e6,1),"M"),n_star = paste0(round(n_star/1e6,1),"M")) %>%
               #knitr::kable("simple") %>% 
               print()
            print("************************************************************")
            
            res %>% mutate(chromosome = strsplit(chromosome,"chr")[[1]][2], alpha = alpha, Age.Only=Age.Only) %>%
               rbind(all_res) -> all_res
               
         }
         
      }
      print("################################################################################################")
   }
   return(all_res)
   
}

get_graph <- function(all_res, alpha_val, dmr=FALSE, scales = "fixed", Age.only=TRUE){
   if(dmr == TRUE) {
      all_res$RMSE[scenario] = all_res$RMSE_dmr
   }
   title = paste0("Alpha = ",alpha_val, "; RMSE by chromosome and model")

   all_res %>% filter(scenario == 3) %>%
      mutate(RMSE = RMSE_dmr, R2 = R2_dmr, scenario = "3 (dmr)") %>%
      rbind(all_res) %>%
      mutate(scenario = ifelse(scenario == 3, "3 (all)", scenario),
             scenario = ifelse(scenario == 4, "4 (dmr)", scenario),
             scenario = ifelse(scenario == 1, "1 (all)", scenario),
             RMSE = ifelse(RMSE > 4000, 0.12, RMSE)) %>%
      mutate(chromosome = as.factor(as.numeric(chromosome))) %>% 
      filter(alpha == alpha_val, Age.Only == Age.only) %>% 
      ggplot(aes(x=chromosome, y=RMSE, color=model, group=model)) +
      geom_point(size=3, shape=19) +  
      geom_line() +  
      facet_wrap(~ scenario, scales = scales,labeller = labeller(scenario = label_both)) +  
      theme_minimal() +  
      theme(
         strip.background = element_blank(),
         strip.text.x = element_text(face = "bold", color = "navy"), 
         legend.position = "bottom",  
         axis.text.x = element_text(face="bold", angle = 0, hjust = 1), 
         panel.spacing = unit(1, "lines") 
      ) +
      labs(
         title = title,
         subtitle = "Faceted by scenario",
         x = "Chromosome",
         y = "RMSE",
         color = "Model"
      ) -> p
   return(p)
}

# run:
setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")
suppressMessages(source("new_code/load_files.R"))
results = get_results_data()
get_graph(results, alpha_val =  .05,dmr=FALSE)
get_graph(results, alpha_val =  .1,dmr=FALSE)
get_graph(results, alpha_val =  .2 ,dmr=FALSE)

