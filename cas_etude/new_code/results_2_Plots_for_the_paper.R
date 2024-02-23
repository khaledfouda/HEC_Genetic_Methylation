setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude")
source("new_code/load_files.R")
library(xtable)

#------------------------------------------
chromosomes = paste0("chr", c(1:12,17))
note = "_subset_Blood"
alphas = c(1e-4)
correction =  function(alpha,N) alpha
#----------------------------------------------
# Table 1
dat.info.full <- read_tsv("./new_data/samples.tsv", show_col_types = FALSE) %>% 
   mutate(bedFile = paste0("new_data/Gene_data/", bedFile)) %>%
   select(DONOR_AGE, bedFile)


set.seed(2022);readRDS("new_data/sample_info.rds") %>% 
   merge(dat.info.full, by='bedFile') %>%
   mutate(DONOR_SEX = replace_na(DONOR_SEX, "Male")) %>% 
   sample_frac(replace = FALSE) %>% 
   group_by(DISEASE, DONOR_SEX) %>%
   filter(row_number() < 5) %>% 
   ungroup() %>% 
   sample_n(15)  %>%
   #dplyr::rename(AGE = DONOR_AGE, SEX=DONOR_SEX) %>% 
   transmute(AGE = DONOR_AGE.y, SEX=DONOR_SEX, DISEASE = DISEASE, TISSUE_TYPE=TISSUE_TYPE,
            CELL_TYPE = mapply(function(x) stringr::str_trunc(x, 31),as.character(CELL_TYPE))) ->
   # select(-bedFile, -featherFile, -XY.ratio, - DONOR_HEALTH_STATUS,
   #        -cellTypeShort,
   #        -sampleGroup) ->
   dat.tab1

tab1 <- xtable(dat.tab1, 
               caption=paste0("Overview of a sample of the donor information with ",
                              "selected columns"),
               label="case2:tab1")
               #align=c("l", "p{3cm}", "p{2cm}", "p{2cm}", "p{3cm}", "p{4cm}", "p{3cm}", "p{3cm}"))
print(tab1, include.rownames = FALSE, caption.placement = "top",
      tabular.environment = 'tabularx', 
      width = "\\textwidth")
#--------------------------------------------------------------------------
# table 2:
# 
# read.csv("new_data/sample_information.csv") %>%
#    select(DONOR_ID, DONOR_AGE) %>% right_join(dat.tab1, by="DONOR_ID",
#                                               relationship="many-to-many")
#--------------------------------------------------

readRDS("new_data/Xdat_common_chr1_subset_Blood.rds") %>% 
   as.data.frame -> x #%>% as.data.frame %>%


x %>% ggplot(aes(x = as.factor(AGE))) +
   geom_histogram(aes(y = ..count..),
                  stat="count", binwidth = 1, fill = "#FFA07A", color = "#FF4500") +  
   theme_minimal() +
   labs(x = "AGE", y = "Count") +
   #scale_fill_manual(values = c("lmcc" = "#FFA07A", "mcc" = "#FFD700")) + 
   theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
   ) +
   ggtitle("Distribution of age among selected samples") -> p1; p1


ggsave("case2_fig1.png", p1, width = 5, height = 4, dpi = 300)
#--------------------------------------------------------------------------
# Fig2: 
dmr.info <- readRDS(paste0("new_data/dmr_info_",note, ".rds"))
dmr.info %>% group_by(chromosome) %>% summarise_all(length) %>% summarise_all(mean)


dmr.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   filter(region_length != 0) %>% 
   group_by(chromosome) %>% 
   summarise(`Number of Regions` = length(site), 
             `Total Number of Sites in the Regions` = sum(region_length)) %>%
   ungroup() %>%
   arrange(chromosome) %>%
   mutate(chromosome = as.factor(chromosome)) %>%
   pivot_longer(cols=c(2,3)) %>%
   ggplot(aes(x = chromosome, y = value, fill = name)) +
   geom_bar(stat = "identity") +  # Draw the bars
   theme_minimal() +  # Use a minimal theme
   scale_fill_manual(values = c("#FFA07A", "#FFD07A")) +  
   facet_wrap(~name, nrow = 2, scales = "free_y") + 
   labs(
      x = "Chromosome", 
      y = "Count",
      title = "Correlated Regions by Chromosome"
      #title = "Distribution of the number of regions and the number of sites within the regions per chromosomes"
   ) +
   theme(
      #axis.text = element_text(face="bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
      axis.title.x = element_text(size = 12, face = "bold"), 
      axis.title.y = element_text(size = 12, face = "bold"), 
      strip.background = element_blank(), 
      strip.text.x = element_text(size = 10, face = "bold"),  
      panel.spacing = unit(1, "lines"),  
      legend.position = "none"  
   ) -> p2; p2

# Save the plot as an image
ggsave("case2_fig2.png", p2, width = 8, height = 6, dpi = 300)
#-------------------------------------------------------------------------------------
library(forcats) 

dmr.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   filter(region_length != 0) %>% 
   group_by(chromosome) %>%
   summarise(
      mean_length = mean(region_length), 
      sd_length = sd(region_length),
      max_length = max(region_length)
   ) %>%
   ungroup() %>%
   arrange(chromosome) %>%
   mutate(chromosome = factor(chromosome)) %>% #summarise_all(mean) 
   ggplot(aes(x = fct_inorder(chromosome), y = mean_length)) +
   geom_bar(stat = "identity", fill = "#FFA07A") + 
   geom_bar(aes(y = max_length), stat = "identity", fill = "#FFA07A", alpha = 0.5) + 
   geom_errorbar(
      aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), 
      width = 0.2, 
      color = "black"
   ) +
   theme_minimal() +  
   theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 1),  
      legend.position = "none"
   ) +
   labs(
      title = "Number of Sites per Region by Chromosome",
      subtitle = "Black lines represent the mean +/- standard deviation;\nLight bars indicate the maximum number of sites per region;\nDark bars represent the average number of sites per region.",
      x = "Chromosome",
      y = "Number of Sites"
   ) -> p3;p3

ggsave("case2_fig3.png", p3, width = 8, height = 6, dpi = 300) 
#------------------------------------------------------------------------------------------

methyl.info <- readRDS(paste0("new_data/methyl_info_",note, ".rds"))



methyl.info %>%
   mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
   arrange(chromosome) %>%
   mutate(dmr = ifelse(dmr == TRUE, "Sites in significant regions", "Sites in non significant regions")) %>% 
   mutate(chromosome = factor(chromosome)) %>%
   ggplot(aes(x = Methylation, group = chromosome)) +
   geom_density(aes( y=after_stat(scaled)),
                fill = "#FFD07A", color = "black", alpha = 0.5) +  
   facet_wrap(~ dmr, nrow = 2) +
   theme_minimal() +
   labs(
      title = "Distribution of Methylation Levels",
      #subtitle = "Different lines for different chromosomes",
      x = "Methylation Level",
      y = "Density"
   ) +
   theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"), 
      axis.text.x = element_text(angle = 45, hjust = 1),  
      panel.spacing = unit(1, "lines"),
      plot.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none" 
   ) -> p4;p4

ggsave("case2_fig4.png", p4, width = 8, height = 6, dpi = 300)  
#------------------------------------------------------------------------
# figure 5

p_values = readRDS(paste0("new_data/p_values_chr2", note, ".rds"))
sites = readRDS(paste0("new_data/sites_common_chr2", note, ".rds"))
site_order = order(sites) 
sites = sites[site_order]
N = length(sites)

dmr_regions = get_dmr_regions(p_values, sites, alpha= 1e-4, min_freq = 1, return_seq = T,
                              middle_point = T,floor_by = 1e3)
data.frame(x = 1:N, Site= sites, NegLogP = - log10(p_values)) %>% 
   mutate(color = ifelse(Site %in% dmr_regions,"DMR","Not DMR")) ->
   manh.dat2

p_values = readRDS(paste0("new_data/p_values_chr3", note, ".rds"))
sites = readRDS(paste0("new_data/sites_common_chr3", note, ".rds"))
site_order = order(sites) 
sites = sites[site_order]
N = length(sites)

dmr_regions = get_dmr_regions(p_values, sites, alpha= 1e-4, min_freq = 1, return_seq = T,
                              middle_point = T,floor_by = 1e3)
data.frame(x = 1:N, Site= sites, NegLogP = - log10(p_values)) %>% 
   mutate(color = ifelse(Site %in% dmr_regions,"DMR","Not DMR")) ->
   manh.dat3



#---------------------------------------------------------------------
#--------------------------------------------------------------------


gtitle <- "Manhattan Plot for Methylation Sites at Chromosome 2"
gnote <- paste("Grey line indicates significance threshold at", expression(1e-4))

manh.dat2 %>% 
   ggplot(aes(x = Site, y = NegLogP, color = color)) +  
   geom_point( size=0.5) +
   scale_color_manual(values = c( "#cc0000","#FFD07A"), 
                      labels = c("Correlated Regions","Uncorrelated Regions")) + 
   theme_minimal() +
   labs(
      x = "Methylation Site",
      y = "- Log (P-value)",#expression("-log"[10] * "(P-value)"),
      title =  gtitle,
      subtitle = gnote
   ) +
   guides(colour = guide_legend(override.aes = list(size=2))) +
   geom_hline(yintercept = -log10(1e-4), color = "gray",alpha=0.4) +
   theme(panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      
      legend.position = "top"  ,
      legend.title =  element_blank(),
      legend.text =  element_text(size = 12, face = "bold")
   ) -> p5; #p5

gtitle <- "Manhattan Plot for Methylation Sites at Chromosome 3"

manh.dat3 %>% 
   ggplot(aes(x = Site, y = NegLogP, color = color)) +  
   geom_point( size=0.5) +
   scale_color_manual(values = c( "#cc0000","#FFD07A"), 
                      labels = c("Correlated Regions","Uncorrelated Regions")) + 
   theme_minimal() +
   labs(
      x = "Methylation Site",
      y = "- Log (P-value)", #expression("-log"[10] * "(P-value)"),
      title =  gtitle,
      subtitle = gnote
   ) +
   guides(colour = guide_legend(override.aes = list(size=2))) +
   geom_hline(yintercept = -log10(1e-4), color = "gray",alpha=0.4) +
   theme(panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
         axis.title.x = element_text(size = 12, face = "bold"),
         axis.title.y = element_text(size = 12, face = "bold"),
         legend.position = "top"  ,
         legend.title =  element_blank(),
         legend.text =  element_text(size = 12, face = "bold")
   ) -> p6; #p6



ggsave("case2_fig5_1.png", p5, width = 8, height = 6, dpi = 300)
ggsave("case2_fig5_2.png", p6, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------------------------


 