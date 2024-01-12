setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")

library(tidyverse)
library(gridExtra)

min_freq = 1

all_res <-data.frame()
   
   
for(Age.Only in c(T,F)){
   for(alpha in c(.05,.1,.2)){
      for(chromosome in paste0("chr",c(7,8,11,12,17))){
         
         if(chromosome == "chr17" & Age.Only == FALSE & alpha == 0.1) next
         print("***********************************************************************")
         load(paste0("results/res_dmr_",chromosome,"_MaleOnly_TRUE_AgeOnly_",Age.Only ,"_",
                     alpha,"x",min_freq,".Rdata"))
         print(paste0("Chromosome ", strsplit(chromosome,"chr")[[1]][2], "; MALE ONLY;",
                      ifelse(Age.Only, " AGE ONLY", " ALL COVARIATES"), "; alpha = ",alpha,"; min_freq = ",min_freq,">>>>>>>>"))
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
}


#res

 # all_res %>%
 #   mutate(chromosome = as.numeric(chromosome)) %>% 
 #   filter(scenario == 4, Age.Only == TRUE) %>% 
 #   ggplot(aes(x=as.factor(chromosome), y=RMSE_dmr, color=model, group=model)) +
 #   geom_point() +
 #   geom_path() +
 #   facet_wrap(~ alpha) 


# library(ggplot2)
# library(dplyr)
# 
# 
# 
# all_res %>%
#    mutate(chromosome = as.factor(as.numeric(chromosome))) %>% 
#    filter(scenario == 4, Age.Only == TRUE) %>% 
#    ggplot(aes(x=chromosome, y=RMSE_dmr, color=model, group=model)) +
#    geom_point(size=3, shape=19) +  
#    geom_line() +  
#    facet_wrap(~ alpha, scales = "fixed") +  
#    theme_minimal() +  
#    theme(
#       strip.background = element_blank(),
#       strip.text.x = element_text(face = "bold", color = "navy"), 
#       legend.position = "bottom",  
#       axis.text.x = element_text(face="bold", angle = 0, hjust = 1), 
#       panel.spacing = unit(1, "lines") 
#    ) +
#    labs(
#       title = "Scenario 4: RMSE on all missing (dmr)",
#       x = "Chromosome",
#       y = "RMSE",
#       color = "Model"
#    ) #-> p1
# 
# 
# all_res %>%
#    mutate(chromosome = as.factor(as.numeric(chromosome))) %>% 
#    filter(scenario == 1, Age.Only == TRUE) %>% 
#    ggplot(aes(x=chromosome, y=RMSE, color=model, group=model)) +
#    geom_point(size=3, shape=19) +  
#    geom_line() +  
#    facet_wrap(~ alpha, scales = "fixed") +  
#    theme_minimal() +  
#    theme(
#       strip.background = element_blank(),
#       strip.text.x = element_text(face = "bold", color = "navy"), 
#       legend.position = "bottom",  
#       axis.text.x = element_text(face="bold", angle = 0, hjust = 1), 
#       panel.spacing = unit(1, "lines") 
#    ) +
#    labs(
#       title = "Scenario 1: RMSE on all missing",
#       x = "Chromosome",
#       y = "RMSE",
#       color = "Model"
#    ) -> p4
# 
# 
# 
# all_res %>%
#    mutate(chromosome = as.factor(as.numeric(chromosome))) %>% 
#    filter(scenario == 3, Age.Only == TRUE) %>% 
#    ggplot(aes(x=chromosome, y=RMSE_dmr, color=model, group=model)) +
#    geom_point(size=3, shape=19) +  
#    geom_line() +  
#    facet_wrap(~ alpha, scales = "fixed") +  
#    theme_minimal() +  
#    theme(
#       strip.background = element_blank(),
#       strip.text.x = element_text(face = "bold", color = "navy"), 
#       legend.position = "bottom",  
#       axis.text.x = element_text(face="bold", angle = 0, hjust = 1), 
#       panel.spacing = unit(1, "lines") 
#    ) +
#    labs(
#       title = "Scenario 3: RMSE on DMR",
#       x = "Chromosome",
#       y = "RMSE (dmr)",
#       color = "Model"
#    ) -> p2
# 
# 
# 
# all_res %>%
#    mutate(chromosome = as.factor(as.numeric(chromosome))) %>% 
#    filter(scenario == 3, Age.Only == TRUE) %>% 
#    ggplot(aes(x=chromosome, y=RMSE, color=model, group=model)) +
#    geom_point(size=3, shape=19) +  
#    geom_line() +  
#    facet_wrap(~ alpha, scales = "free_y") +  
#    theme_minimal() +  
#    theme(
#       strip.background = element_blank(),
#       strip.text.x = element_text(face = "bold", color = "navy"), 
#       legend.position = "bottom",  
#       axis.text.x = element_text(face="bold", angle = 0, hjust = 1), 
#       panel.spacing = unit(1, "lines") 
#    ) +
#    labs(
#       title = "Scenario 3: RMSE on all missing",
#       x = "Chromosome",
#       y = "RMSE",
#       color = "Model"
#    ) -> p3


all_res %>% filter(scenario == 3) %>%
   mutate(RMSE = RMSE_dmr, R2 = R2_dmr, scenario = "3 (dmr)") %>%
   rbind(all_res) %>%
   mutate(scenario = ifelse(scenario == 3, "3 (all)", scenario),
          scenario = ifelse(scenario == 4, "4 (dmr)", scenario),
          scenario = ifelse(scenario == 1, "1 (all)", scenario),
          RMSE = ifelse(RMSE > 4000, 0.12, RMSE)) ->
   new_res




get_graph <- function(all_res, alpha_val, dmr=TRUE, scales = "fixed"){
   if(dmr == TRUE) all_res$RMSE = all_res$RMSE_dmr
   title = paste0("Alpha = ",alpha_val, "; RMSE by chromosome and model")
   
   all_res %>%
      mutate(chromosome = as.factor(as.numeric(chromosome))) %>% 
      filter(alpha == alpha_val, Age.Only == TRUE) %>% 
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

get_graph(new_res, .05, F, "fixed")
get_graph(new_res, .1, F, "fixed")
get_graph(new_res, .2, F, "fixed")

# grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# 
# all_res %>%
#    mutate(chromosome = as.numeric(chromosome)) %>% 
#    filter(scenario == 3, Age.Only == TRUE) %>% 
#    ggplot(aes(x=as.factor(chromosome), y=RMSE_dmr, color=model, group=model)) +
#    geom_point() +
#    geom_path() +
#    facet_wrap(~ alpha)
# 
# 
# all_res %>%
#    mutate(chromosome = as.numeric(chromosome)) %>% 
#    filter(scenario == 3, Age.Only == TRUE) %>% 
#    ggplot(aes(x=as.factor(chromosome), y=RMSE, color=model, group=model)) +
#    geom_point() +
#    geom_path() +
#    facet_wrap(~ alpha, scales = "free_y")
# 
# 
# all_res %>%
#    mutate(chromosome = as.numeric(chromosome)) %>% 
#    filter(scenario == 1, Age.Only == TRUE) %>% 
#    ggplot(aes(x=as.factor(chromosome), y=RMSE, color=model, group=model)) +
#    geom_point() +
#    geom_path() +
#    facet_wrap(~ alpha)
# 
# 
# p <- ggplot(res, aes(x = factor(scenario), y = RMSE, group = model, color=model, fill=model)) + geom_point()
# ggplotly(p)
# p
# 
# fig <- plot_ly(data, x = ~x, y = ~y, type = 'scatter', mode = 'markers')
# fig
