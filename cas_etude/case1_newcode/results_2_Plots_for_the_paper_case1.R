setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude")
source("new_code/load_files.R")
library(xtable)

#------------------------------------------
alphas = c(1e-10)
correction =  function(alpha, N)
   alpha
note = "alphae10"
#----------------------------------------------

#-------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Fig2:
#\
dmr.info <-
   readRDS(paste0("new_data/case1_dmr_info_", note, ".rds"))
# need to recreate this
# <site, region_length>



dmr.info %>% summarise_all(length) %>% summarise_all(mean)

head(dmr.info)

dmr.info %>%
   filter(region_length != 0) %>%
   summarise(
      `Number of Regions` = length(site),
      `Total Number of Sites in the Regions` = sum(region_length)
   ) %>%
   mutate(index = 1) %>%
   pivot_longer(cols = c(1, 2)) #%>%
   ggplot(aes(x = index, y = value, fill = name)) +
   geom_bar(stat = "identity") +  # Draw the bars
   theme_minimal() +  # Use a minimal theme
   scale_fill_manual(values = c("#FFA07A", "#FFD07A")) +
   facet_wrap( ~ name, nrow = 2, scales = "free_y") +
   labs(x = "Chromosome",
        y = "Count",
        title = "Correlated Regions by Chromosome") +
   #title = "Distribution of the number of regions and the number of sites within the regions per chromosomes") +
   theme(
      #axis.text = element_text(face="bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10, face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
   ) -> p2
p2

# Save the plot as an image
ggsave(
   "case1_fig2.png",
   p2,
   width = 8,
   height = 6,
   dpi = 300
)
#-------------------------------------------------------------------------------------
library(forcats)

dmr.info %>%
   mutate(index = 1) %>%
   filter(region_length != 0) %>%
   group_by(index) %>%
   summarise(
      mean_length = mean(region_length),
      sd_length = sd(region_length),
      max_length = max(region_length)
   ) %>%
   ungroup() #%>%
   
   ggplot(aes(x = index, y = mean_length)) +
   geom_bar(stat = "identity", fill = "#FFA07A") +
   geom_bar(
      aes(y = max_length),
      stat = "identity",
      fill = "#FFA07A",
      alpha = 0.5
   ) +
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
   ) -> p3
p3

ggsave(
   "case1_fig3.png",
   p3,
   width = 8,
   height = 6,
   dpi = 300
)
#------------------------------------------------------------------------------------------

methyl.info <-
   readRDS(paste0("new_data/case1_methyl_info_", note, ".rds"))



methyl.info %>%
   mutate(index = 1) %>%
   mutate(dmr = ifelse(
      dmr == TRUE,
      "DMR sites", #"Sites in significant regions",
      "non-DMR sites"#"Sites in non significant regions"
   )) %>%
   ggplot(aes(x = Methylation)) +
   geom_density(
      aes(y = after_stat(scaled)),
      fill = "#FFD07A",
      color = "black",
      alpha = 0.5
   ) +
   facet_wrap( ~ dmr, nrow = 2) +
   theme_minimal() +
   labs(title = "Distribution of Methylation Levels",
        #subtitle = "Different lines for different chromosomes",
        x = "Methylation Level",
        y = "Density") +
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

ggsave(
   "case1_fig4.png",
   p4,
   width = 8,
   height = 6,
   dpi = 300
)
#------------------------------------------------------------------------
# figure 5

site_order = order(sites)
sites = sites[site_order]
N = length(sites)

data.frame(x = 1:N,
           Site = sites,
           NegLogP = -log10(p_values)) %>%
   mutate(color = ifelse(Site %in% dmr_regions, "DMR", "Not DMR")) ->
   manh.dat2

#---------------------------------------------------------------------
#--------------------------------------------------------------------


gtitle <-
   "Manhattan Plot for Methylation Levels"
gnote <-
   paste("Grey line indicates significance threshold at",
         expression(1e-10))

manh.dat2 %>%
   ggplot(aes(x = Site, y = NegLogP, color = color)) +
   geom_point(size = 0.5) +
   scale_color_manual(
      values = c("#cc0000", "#FFD07A"),
      labels = c("DMR Regions", "non-DMR Regions")
   ) +
   theme_minimal() +
   labs(
      x = "Methylation Site",
      y = "- Log (P-value)",
      #expression("-log"[10] * "(P-value)"),
      title =  gtitle,
      subtitle = gnote
   ) +
   guides(colour = guide_legend(override.aes = list(size = 2))) +
   geom_hline(
      yintercept = -log10(alpha),
      color = "gray",
      alpha = 0.4
   ) +
   theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      
      legend.position = "top"  ,
      legend.title =  element_blank(),
      legend.text =  element_text(size = 12, face = "bold")
   ) -> p5;p5




ggsave(
   "case1_fig5_1.png",
   p5,
   width = 8,
   height = 6,
   dpi = 300
)

#-------------------------------------------------------------------------------------------------
# Load necessary package
library(xtable)

load(paste0(
   "results/res_case1_",
   alpha, note, ".Rdata"))
res

res %>% 
   mutate(scenario = paste("scenario", scenario),
          Computation_time = time / 60) %>%
   dplyr::rename(Method = model,
          Q2 = R2,
          Q2_dmr = R2_dmr) %>%
   mutate_if(is.numeric, round, digits=3) %>%  
   dplyr::select(scenario, Method, RMSE, Q2, RMSE_dmr,
          Q2_dmr, Computation_time) -> data

head(data)

data <- data.frame(
   Scenario = rep(c(
      "Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"
   ), each = 3),
   Method = rep(c("LMCC", "GaSP", "Naive"), 4),
   RMSE = runif(12, 0.1, 0.3),
   Q2 = runif(12, 0.1, 0.9),
   RMSE_dmr = c(rep(NA, 3), runif(9, 0.1, 0.3)),
   Q2_dmr = c(rep(NA, 3), runif(9, 0.1, 0.9)),
   Computation_time = runif(12, 0.1, 35)
)

# Convert NAs to '-' for LaTeX compatibility
data[is.na(data)] <- "-"

# Create xtable object from the dataframe
xt <- xtable(data,
             caption = paste("Prediction results (RMSE and Q2) and computational time", 
                             "for the methylation data set. \\textbf{Scenario 1}: NA's",
                             "are uniformly distributed on non DMR sites (90,000 sites).",
                             "\\textbf{Scenario 2}: NA's are partially located on some",
                             "potential DMR's (2,000 sites) and the remaining sites are",
                             "located on non DMR sites (90,000 sites).",
                             "\\textbf{Scenario 3}: NA's are located on all potential",
                             "DMR's (4,000 sites) and the remaining sites are located on",
                             "non DMR sites (90,000 sites). \\textbf{Scenario 4}: NA's",
                             "are located only on potential DMR's (4000 sites)"),
             label = "tab:methyl2")

# Print the LaTeX table code
print(
   xt,
   include.rownames = FALSE,
   hline.after = c(-1, 0, 3, 6, 9),#12
   sanitize.text.function = function(x) {
      x
   }
)
