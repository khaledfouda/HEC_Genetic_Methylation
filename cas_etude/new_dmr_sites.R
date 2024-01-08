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


Y1 = readRDS("new_data/Ydat_common_chr1.rds")
X1 = readRDS("new_data/Xdat_common_chr1.rds") %>% as.data.frame() #%>% 
sites1 = readRDS("new_data/sites_common_chr1.rds")

Y2 = readRDS("new_data/Ydat_common_chr2.rds")
X2 = readRDS("new_data/Xdat_common_chr2.rds") %>% as.data.frame() #%>% 
sites2 = readRDS("new_data/sites_common_chr2.rds")
all(X1[,1:4]==X2[,1:4])
Y = cbind(Y1,Y2)
sites = c(sites1, sites2)
sites_order = order(sites)
Y = Y[,sites_order]
sites = sites[sites_order]
X=X1
dim(Y)
   #mutate(AGE = (AGE - mean(AGE))/ sd(AGE) ) %>% 
   #as.matrix()
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
   xlab("Methylation Site") +
   ylab("-log10(Adjusted P-value)") +
   geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red")+
   #scale_x_continuous(breaks = manh.dat$x, labels = manh.dat$sites)  +
   ggtitle("Manhattan Plot of Methylation Sites") -> p 
   
ggsave(filename = "manhattan_plot_chr1_chr2.png", plot = p, width = 10, height = 6, dpi = 300)
p
 
range(p_values,na.rm = T) 


wherenais = which(is.na(p_values))
length(wherenais)
length(p_values)
#-------------------------------------------------------

model <- lm(Y[, 1] ~ X$AGE)
summary(model)$coefficients["X$AGE", "Pr(>|t|)"]


p.adjust(c(.1,.05,.3), method="bonferroni", n=3)

N

n_star.new = sort(unique(c(round(seq(1,Nnew,length.out=round(Nnew*0.90))))))
sites.new = scale_01(1:Nnew)

k_star.new <-
   (Xnew %>%
       as.data.frame() %>%
       mutate(index = 1:n()) %>%  
       group_by(MALE, AML, APL, BONE_MARROW) %>% 
       filter(row_number() != 1) %>% 
       ungroup())$index  




