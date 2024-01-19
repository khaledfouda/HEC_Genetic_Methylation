

setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")  



library(tidyverse)

chromosome = "chr2"
note = "_subset_Blood"
correction <- function(alpha,N) alpha
alpha = 1e-4
region_length <- 1e3
dmr.info <- methyl.info <- data.frame()

run_code = FALSE
if(run_code){
   for(chromosome in paste0("chr",c(1:12,17))){
      print(chromosome)
      sites = readRDS(paste0("new_data/sites_common_",chromosome, note, ".rds"))
      p_values = readRDS(paste0("new_data/p_values_",chromosome, note, ".rds"))
      N = length(sites)
      alpha = correction(alpha, N)
      loci = sort(sites[which(p_values < alpha)])
      sequen = numeric(length(loci))
      step_size = round(region_length/2)
      end_site = 1
      
      for(i in 1:length(loci)){
         start_site = max( end_site, loci[i]-step_size)
         end_site = loci[i] + step_size
         sequen[i] = length(which(sites %in%  (start_site: end_site)))
      }
      
      data.frame(chromosome = chromosome, site=loci, region_length=sequen) %>%
         rbind(dmr.info) ->
         dmr.info
      
      #-------------------------
      # Get Methylation info
      dmr_regions = get_dmr_regions(p_values, sites, alpha, floor_by = region_length,
                                     middle_point = T, return_seq = T)
      ind_dmr = sort(which(sites %in% dmr_regions))
      #sites = sites[ind_dmr]
      Y = readRDS(paste0("new_data/Ydat_common_",chromosome, note, ".rds"))
      Y = colMeans(Y)
      # data.frame(chromosome = chromosome, site = sites, Methylation = Y[ind_dmr]) %>%
      data.frame(chromosome = chromosome, site = sites, Methylation = Y, dmr = 1:N) %>%
         mutate(dmr = ifelse(dmr %in% ind_dmr, TRUE, FALSE)) %>% 
         rbind(methyl.info) ->
         methyl.info
      
   }
   saveRDS(dmr.info, paste0("new_data/dmr_info_",note, ".rds"))
   saveRDS(methyl.info, paste0("new_data/methyl_info_",note, ".rds"))
   
}
dmr.info <- readRDS(paste0("new_data/dmr_info_",note, ".rds"))
methyl.info <- readRDS(paste0("new_data/methyl_info_",note, ".rds"))


dmr.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   filter(region_length != 0) %>% 
   group_by(chromosome) %>% 
   summarise(num_DMR = length(site), 
             num_sites_in_DMR = sum(region_length)) %>%
   ungroup() %>%
   arrange(chromosome) %>%
   mutate(chromosome = as.factor(chromosome)) %>%
   kable(format = "pipe")


methyl.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   arrange(chromosome) %>%
   mutate(chromosome = as.factor(chromosome)) %>%
   #filter(chromosome == 17) %>% 
   ggplot(aes(group=chromosome)) +
   geom_histogram(aes(Methylation,y=after_stat(scaled))) +
   facet_wrap(~ dmr, nrow=2)



methyl.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   arrange(chromosome) %>%
   #sample_n(1e6) %>% 
   mutate(dmr = ifelse(dmr==TRUE, "DMR", "Not DMR")) %>% 
   mutate(chromosome = as.factor(chromosome)) %>%
   ggplot(aes(x = Methylation, group = chromosome)) +
   #geom_histogram( bins = 30, fill = "steelblue", color = "black") +
   geom_density(aes(y=after_stat(scaled)), fill="steelblue", color="black",alpha=0.01) +
   facet_wrap(~ dmr, nrow = 2) +
   theme_minimal() +
   labs(
      title = "Methylation Distribution Comparison between DMR and non-DMR sites",
      subtitle = "Different lines for different chromosomes",
      x = "Methylation Level",
      y = "Density"
   ) +
   theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"), 
      axis.text.x = element_text(angle = 0, hjust = 1),
      panel.spacing = unit(1, "lines"),
      plot.background = element_blank(),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
   )


library(dplyr)
library(ggplot2)

methyl.info %>%
   filter(dmr == TRUE) %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   arrange(desc(chromosome)) %>%
   mutate(chromosome = as.factor(chromosome)) %>%
   ggplot(aes(x = site, y = chromosome)) +
   geom_point(alpha = 1, color = "steelblue", size = 1) +  
   theme_minimal() +
   labs(
      title = "DMRs' Sites location per Chromosome",
      subtitle = "To visually identify common DMR sites among chromosomes, if any.",
      x = "Site",
      y = "Chromosome"
   ) +
   theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text.y = element_text(angle = 0, vjust = 1),
      panel.grid = element_blank()
   )


methyl.info %>%
   filter(dmr == TRUE) %>% 
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   arrange(chromosome) %>%
   mutate(chromosome = as.factor(chromosome)) %>%
   ggplot(aes(x = site, y= chromosome)) +
   geom_point()
   #geom_histogram( bins = 30, fill = "steelblue", color = "black") 
   #geom_density(aes(y=after_stat(scaled)), fill="steelblue", color="black",alpha=0.01) +
   

plot(density(Y))


library(dplyr)
library(ggplot2)

dmr.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   filter(region_length != 0, chromosome != 9) %>%
   ggplot(aes(x = region_length)) +
   geom_histogram(bins = 100, fill = "steelblue", color = "black") +
   theme_minimal() +
   labs(
      title = "Distribution of Number of Sites per Single DMR",
      x = "Number of Sites",
      y = "Count"
   ) +
   theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, colour = "steelblue", face="bold"),
      axis.title.y = element_text(size = 12, colour = "steelblue", face="bold"),
      axis.text.x = element_text(angle = 0, hjust = 1, face = "bold", color="red"),
      axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", color="red")
   ) +
   facet_wrap(~chromosome, labeller = label_both)




dmr.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   filter(region_length != 0) %>% 
   group_by(chromosome) %>%
   summarise(mean_length = mean(region_length), 
             sd_length = sd(region_length),
             max_length = max(region_length)) %>%
   ungroup() %>%
   arrange(chromosome) %>%
   mutate(chromosome = as.factor(chromosome)) %>% 
   ggplot(aes(x = fct_inorder(chromosome), y = mean_length)) +
   geom_bar(stat = "identity", fill = "steelblue") +
   geom_bar(aes(y=max_length), stat = "identity", fill = "steelblue", alpha=0.5) +
   geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), width = 0.2) +
   theme_bw() +
   theme(
      axis.text.x = element_text(angle = 0, hjust = 1), 
      axis.title = element_text(size = 12),
      title = element_text(size = 14)
   ) +
   labs(
      title = "Region Length by Chromosome",
      subtitle = "Black lines are the mean +/- standard deviation;\nLight bars are the Max lengths;\nDark bars the average lengths.",
      x = "Chromosome",
      y = "Mean / Max Region Length"
   )


dmr.info %>%
   filter(region_length != 0) %>% 
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   group_by(chromosome) %>%
   summarise(num_p_values = length(region_length),
             num_sites = sum(region_length)) %>%
   pivot_longer(cols = -chromosome, names_to = "statistic", values_to = "Value") %>%
   mutate(statistic = ifelse(statistic == "num_sites", "Number of Sites","Number of p-vales")) %>% 
   arrange(chromosome) %>% 
   mutate(chromosome = as.factor(chromosome)) %>% 
   ggplot(aes(x = fct_inorder(chromosome), y = Value)) +
   geom_bar(stat = "identity", fill = "steelblue") +
   facet_wrap(~statistic, nrow=2, scales = "free_y") +   
   theme_bw() +
   theme(
      axis.text.x = element_text(angle = 0, hjust = 1), 
      axis.title = element_text(size = 12),
      title = element_text(size = 14)
   ) +
   labs(
      title = "Number of sites/p-values for each Chromosome",
      subtitle = paste0("Number of sites: total number of sites chosen as DMR;",
                        "\nNumber of p-values: Total number of p-values exceeding the threshold."),
      x = "Chromosome",
      y = ""
   )

dmr.info %>%
   filter(region_length != 0) %>% 
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   #filter(chromosome == 2) %>% 
   arrange(chromosome) %>% 
   mutate(chromosome = as.factor(chromosome)) %>% 
   ggplot(aes(fill=chromosome)) +
   geom_density(aes(x=region_length, y=after_stat(scaled)), alpha=0.3) +
   theme_bw() +
   theme(
      axis.text.x = element_text(angle = 0, hjust = 1), 
      axis.title = element_text(size = 12),
      title = element_text(size = 14)
   ) +
   labs(
      title = "Distribution of Region Length Grouped by Chromosomes",
      subtitle = paste0("The region of sites per chosen p-value."),
      x = "Region Length",
      y = "Density"
   )


# 
# dmr.info %>% 
#    filter(region_length != 0) %>% 
#    group_by(chromosome) %>%
#    summarise(
#       Min = min(region_length),
#       Median = median(region_length),
#       Mean = mean(region_length),
#       Max = max(region_length),
#       SD = sd(region_length)
#    )
# 
# dmr.info %>%
#    filter(region_length == 0)
# 
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# 
# stats <- dmr.info %>%
#    filter(region_length != 0) %>% 
#    mutate(chromosome = factor(gsub("chr", "", chromosome), levels = unique(gsub("chr", "", chromosome)))) %>%
#    group_by(chromosome) %>%
#    summarise(
#       Average = mean(region_length),
#       Max = max(region_length),
#       Median = median(region_length)
#    ) %>%
#    ungroup()
# 
# long_stats <- stats %>%
#    pivot_longer(cols = -chromosome, names_to = "Statistic", values_to = "Value")
# 
# ggplot(long_stats, aes(x = chromosome, y = Value, group = Statistic, color = Statistic)) +
#    geom_line() +
#    theme_minimal() +
#    theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       axis.title = element_text(size = 12),
#       title = element_text(size = 14)
#    ) +
#    labs(
#       title = "Chromosomal Statistics of Region Length",
#       x = "Chromosome",
#       y = "Region Length",
#       color = "Statistic"
#    )
# 
# 




