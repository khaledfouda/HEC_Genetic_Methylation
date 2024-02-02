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


set.seed(2022);readRDS("new_data/sample_info.rds") %>%
   mutate(DONOR_SEX = replace_na(DONOR_SEX, "Male")) %>% 
   sample_frac(replace = FALSE) %>% 
   group_by(DISEASE, DONOR_SEX) %>%
   filter(row_number() < 5) %>% 
   ungroup() %>% 
   sample_n(15)  %>%
   #dplyr::rename(AGE = DONOR_AGE, SEX=DONOR_SEX) %>% 
   mutate(CELL_TYPE = mapply(function(x) stringr::str_trunc(x, 31),as.character(CELL_TYPE))) %>% 
   select(-bedFile, -featherFile, -XY.ratio, - DONOR_HEALTH_STATUS,
          -cellTypeShort,
          -sampleGroup) ->
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
methyl.info <- readRDS(paste0("new_data/methyl_info_",note, ".rds"))


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
   scale_fill_manual(values = c("#FFA07A", "#FFD700")) +  
   facet_wrap(~name, nrow = 2, scales = "free_y") + 
   labs(
      x = "Chromosome", 
      y = "Count",
      title = "Distribution of the number of regions and the number of sites within the regions per chromosomes"
   ) +
   theme(
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





