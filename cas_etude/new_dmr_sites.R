setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")

library(Jmisc)
library(DiceEval)
library(parallel)
library(plotly)
library(mgcv)
library(Rcpp)
library(RColorBrewer)
library(ggpubr)
library(ZIprop)
library(corrplot)


source("methyl_func.R")
sourceAll(path="../functions/")

chromosome = "chr3"


#---------------
# combined
chromosome = "chr8"
run_p_values_plot("chr12")
run_p_values_plot("chr17")
run_p_values_plot <- function(chromosome){
Y = readRDS(paste0("new_data/Ydat_common_",chromosome,".rds"))
X = readRDS(paste0("new_data/Xdat_common_",chromosome,".rds")) %>% as.data.frame() #%>% 
sites = readRDS(paste0("new_data/sites_common_",chromosome,".rds"))

# for(chromosome in paste0("chr",c(8,11,12,17))){
#    print(chromosome)
#    Y = cbind(Y, readRDS(paste0("new_data/Ydat_common_",chromosome,".rds")))
#    #X = readRDS(paste0("new_data/Xdat_common_",chromosome,".rds")) %>% as.data.frame() #%>% 
#    sites = c(sites,readRDS(paste0("new_data/sites_common_",chromosome,".rds")))
# }

female_indices = which(X$MALE == 1)
Y = Y[female_indices,]
X = X[female_indices,]

site_order = order(sites) 
Y = Y[,site_order]
sites = sites[site_order]

dim(Y)
N = ncol(Y)
K = nrow(Y)
#--------------------------------
library(parallel)
no_cores <- 12               #detectCores() - 1  # Leave one core free for system processes
cl <- makeCluster(no_cores)

Age = X$AGE
clusterExport(cl, varlist = c("Y", "Age"))

p_values <- parSapply(cl, 1:ncol(Y), function(i) {
   model <- lm(Y[, i] ~ Age) 
   summary(model)$coefficients["Age", "Pr(>|t|)"]
}) 
stopCluster(cl)
saveRDS(p_values, paste0("new_data/p_values_",chromosome,".rds"))
#adjusted_p_values <- p.adjust(p_values, method = "bonferroni", n = length(p_values))  
#-------------------------------------------------------------
significance_threshold = -log10(0.1 / N) 
data.frame(x = 1:N, Site= sites, NegLogP = - log10(p_values)) ->
   manh.dat
#---------------------------------------------------------------------
signf.p <- which(p_values < (.1/N))
#length(p_values)
#sites[signf.p]
#sites[which(p_values < (0.2/N))]

get_dmr_regions <- function(p_values, N, alpha=0.1, floor_by=1e7, min_freq=2, return_seq=FALSE){
   as.data.frame(table(floor(sites[which(p_values < (alpha/N))] / floor_by) * floor_by)) %>%
      arrange(desc(Freq)) %>%
      filter(Freq >= min_freq) -> results
   if(return_seq == TRUE){
      values = as.numeric(as.character(results$Var1)) 
      sequen = c()
      for(x in values) sequen = c(sequen, x:(x+1e7))
      return(sequen)
   }
   return(results)
}

#get_dmr_regions(p_values, N, 0.05, min_freq = 2, return_seq = T)
#get_dmr_regions(p_values, N, 0.1)
#get_dmr_regions(p_values, N, 0.2)
dmr_regions = get_dmr_regions(p_values, N, 0.05, min_freq = 1, return_seq = T)
manh.dat %<>% mutate(color = ifelse(Site %in% dmr_regions,"DMR","Not DMR"))
#--------------------------------------------------------------------

manh.dat %>% 
   ggplot(aes(Site, NegLogP, color=color, fill=color)) +  
   geom_point(alpha = 0.6) +
   scale_color_manual(values = c("DMR" = "blue", "Not DMR" = "black")) +
   #theme_minimal() +
   #scale_x_continuous(limits=c(0,2e8))+
   xlab("Methylation Site") +
   ylab("-log10(Adjusted P-value)") +
   geom_hline(yintercept = -log10(0.05/N), linetype = "dashed", color = "red")+ 
   geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red")+ 
   geom_hline(yintercept = -log10(0.2/N), linetype = "dashed", color = "red")+ 
   #scale_x_continuous(breaks = manh.dat$x, labels = manh.dat$sites)  +
   ggtitle(paste0("Chromosomes ",chromosome,"; Male","; Alpha=0.5")) -> p 

ggsave(filename = paste0("./graphs/manhattan_plot_",chromosome,"_Male", ".png"),
       plot = p, width = 10, height = 6, dpi = 300) 
}
#------------------------
for(chromosome in paste0("chr",c(7,8,11,12,17))){
   print(chromosome)
   Y = readRDS(paste0("new_data/Ydat_common_",chromosome,".rds"))
   X = readRDS(paste0("new_data/Xdat_common_",chromosome,".rds")) %>% as.data.frame() #%>% 
   sites = readRDS(paste0("new_data/sites_common_",chromosome,".rds"))
   
   dim(Y)
   N = ncol(Y)
   K = nrow(Y)
   #--------------------------------
   library(parallel)
   no_cores <- 15                #detectCores() - 1  # Leave one core free for system processes
   cl <- makeCluster(no_cores)
   
   Age = X$AGE
   clusterExport(cl, varlist = c("Y", "Age"))
   
   p_values <- parSapply(cl, 1:ncol(Y), function(i) {
      model <- lm(Y[, i] ~ Age)
      summary(model)$coefficients["Age", "Pr(>|t|)"]
   })
   stopCluster(cl)
   adjusted_p_values <- p.adjust(p_values, method = "bonferroni", n = length(p_values))
   #saveRDS(p_values, "new_data/p_values.rds")
   #saveRDS(adjusted_p_values, "new_data/adjusted_p_values.rds")
   #-------------------------------------------------------------
   significance_threshold = -log10(0.1 / N) 
      #mutate(Category = ifelse(log10(p_values) < significance_threshold, "Significant", "Not Significant"))
   
   data.frame(x = 1:N, Site= sites, NegLogP = - log10(p_values)) ->
      manh.dat
   manh.dat %>%  
      ggplot(aes(Site, NegLogP)) +  
      geom_point(alpha = 0.6) +
      #scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
      theme_minimal() +
      scale_x_continuous(limits=c(0,2e8))+
      xlab("Methylation Site") +
      ylab("-log10(Adjusted P-value)") +
      geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red")+
      #scale_x_continuous(breaks = manh.dat$x, labels = manh.dat$sites)  +
      ggtitle(paste0("Chromosome ",strsplit(chromosome,"chr")[[1]][2], ", Alpha=0.1")) -> p 
      
   ggsave(filename = paste0("./graphs/manhattan_plot_",chromosome, ".png"),
          plot = p, width = 10, height = 6, dpi = 300)
}


# p
#  
# range(p_values,na.rm = T)  
# 
# 
# wherenais = which(is.na(p_values))
# length(wherenais)
# length(p_values)
# #-------------------------------------------------------
# 
# model <- lm(Y[, 1] ~ X$AGE)
# summary(model)$coefficients["X$AGE", "Pr(>|t|)"]
# 
# 
# p.adjust(c(.1,.05,.3), method="bonferroni", n=3)
# 
# N
# 
# n_star.new = sort(unique(c(round(seq(1,Nnew,length.out=round(Nnew*0.90))))))
# sites.new = scale_01(1:Nnew)
# 
# k_star.new <-
#    (Xnew %>%
#        as.data.frame() %>%
#        mutate(index = 1:n()) %>%  
#        group_by(MALE, AML, APL, BONE_MARROW) %>% 
#        filter(row_number() != 1) %>% 
#        ungroup())$index  
# 
# 
# 
# 
