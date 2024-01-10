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
library(tidyverse)
library(magrittr)

library(doParallel)
library(foreach)

source("methyl_func.R")
sourceAll(path="../functions/")

# ## get data from fastgasp article
# load("data/Methylation_level_dense_Chr1.Rda")
# ## get identified DMR area
# load(file = "data/area_diff.Rdata")
# load(file = "data/area_hypo.Rdata")
# 
# 
# 
# ind = 1:1e6
# Y = t(seq.data.chr1[ind, ])
# selec = c(1:2, 4:9, 12:14, 16:19)
# 
# Y = Y[selec, ]
# colnames(Y) = rownames(Y) = NULL
# X1 = c(rep("Gifford", 8), rep("roadmap1", 7))
# 
# sites = scale_01(Index_seq[ind])
# 
# 
# N = ncol(Y)
# K = nrow(Y)
# methyl = Y
# k_star = (1:K)[-c(1,9)]
# #Y_mean = apply(Y,2,mean)
# X = fact2mat(X1)

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


#### set n_star ####
# 
# area_all = c(area_diff,area_hypo)
# 
# ### 50 % of DMR area are taken
# dmr_50_random = sort(unlist(lapply(area_all,function(x){
#   sample(x,round(length(x)*0.5))
# })))
# dmr_50_random = which(ind %in% dmr_50_random )
# 
# ### all DMR
# ind_dmr = sort(unique(unlist(area_all)))
# ind_dmr = which(ind %in% ind_dmr )
# 
# 
# 
# n_star_list = list(sort(unique(c(round(seq(1,N,length.out=round(N*0.9)))))), #scenario 1
#                    sort(unique(c(round(seq(1,N,length.out=round(N*0.9))),dmr_50_random))), #scenario 2
#                    sort(unique(c(round(seq(1,N,length.out=round(N*0.9))),ind_dmr))), #scenario 3
#                    ind_dmr) #scenario 4
# 
# save(n_star, file = "data/n_star.Rdata")
# 
# sapply(n_star_list,length)/N
# 
# 
# #### apply method ####
# res = list()
# length(res) = length(n_star_list)
# 
# #i=1
# for ( i in 1:4){
#   n_star = n_star_list[[i]]
#   if( i == 2){
#     ind_na_sub = dmr_50_random
#   }else if ( i %in% c(3,4)){
#     ind_na_sub = ind_dmr
#   }else if ( i == 1){
#     ind_na_sub = NULL
#   }
#   res[[i]] = methyl_func(methyl, sites, k_star, n_star, X,ind_na_sub)
#   save(res, file = "data/res_dmr_1e6_bis.Rdata")
# }  
# 
# res
 

chromosome="chr12"; Age.Only=TRUE; Male.Only=TRUE;alpha=.05; min_freq=1; min_k=2;subset=1e5
#-----------------------
run_model_on_chr <- function(chromosome, Age.Only=TRUE, Male.Only=TRUE,alpha=.05, min_freq=1, min_k=2,subset=NA){
    
  setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")  
  Y = readRDS(paste0("new_data/Ydat_common_",chromosome,".rds"))
  sites = readRDS(paste0("new_data/sites_common_",chromosome,".rds"))
  X = readRDS(paste0("new_data/Xdat_common_",chromosome,".rds")) %>%
    as.data.frame() %>% 
    mutate(AGE = (AGE - mean(AGE))/ sd(AGE) ) %>% 
    as.matrix()
  p_values = readRDS(paste0("new_data/p_values_",chromosome, ".rds"))
  
  if(!is.na(subset)){
    Y = Y[,1:subset]
    sites = sites[1:subset]
    p_values = p_values[1:subset]
  }
  
  if(Male.Only == TRUE){
    male_indices = which(X[,"MALE"]==1)
    Y = Y[male_indices,]
    X = X[male_indices,]
  }
  
  N = ncol(Y)
  K = nrow(Y)
  dmr_regions = get_dmr_regions(p_values, N, alpha, min_freq = min_freq, return_seq = T)
  ind_dmr = sort(which(sites %in% dmr_regions))
  scaled_sites = scale_01(sites)
  
  #n_star = sort(unique(c(round(seq(1,Nnew,length.out=round(Nnew*0.90))))))
  #n_star.new = sort(c(n_star.new, ind_na_sub))
  
  
  k_star <-
    (X %>%
    as.data.frame() %>%
    mutate(index = 1:n()) %>%  
    group_by(MALE, AML, APL, BONE_MARROW) %>% 
    filter(row_number() != 1) %>% 
    ungroup())$index  
  
  if( (K-length(k_star)) < min_k ){
    extra_k = sample(1:length(k_star.new), (min_k-K+length(k_star)),replace = FALSE )
    k_star <- k_star[-extra_k]
  }
  if(Age.Only == TRUE){
    X = X[,"AGE"]
  }
  #Xnew <- as.matrix( sample(0:1, Knew, replace=TRUE)) 
  
  #methyl_func(Ynew, sites.new, k_star.new, n_star.new, Xnew, ind_na_sub)
  
  #---------------
  n_star_list = list(sort(unique(c(round(seq(1,N,length.out=round(N*0.9)))))), #scenario 1,
                     sort(unique(c(round(seq(1,N,length.out=round(N*0.9))),ind_dmr))), #scenario 3
                     ind_dmr) #scenario 4
  
  
  
  
  # Detect the number of cores
  no_cores <- length(n_star_list) #detectCores() - 1
  cl <- makeCluster(no_cores)
  nothing = clusterEvalQ(cl, {
    library(Jmisc)
    source("methyl_func.R")
    sourceAll(path="../functions/")
  })
  clusterExport(cl, c("methyl_func", "Y", "scaled_sites", "k_star", "X", "ind_dmr", "n_star_list"))
  
  res <- parLapply(cl, 1:3, function(i) {
    n_star <- n_star_list[[i]]
    ind_na_sub <- ifelse(i == 1, NULL, ind_dmr)
    methyl_func(Y, scaled_sites, k_star, n_star, X, ind_na_sub)
  })
  stopCluster(cl)
  res$n_star = length(n_star)
  res$k_star = length(k_star)
  res$N = N
  res$K = K
  save(res, file = paste0("results/res_dmr_",chromosome,"_MaleOnly_",Male.Only,"_AgeOnly_",
                          Age.Only,"_",alpha,"x",min_freq,".Rdata"))
  return(res)
}

scaled_sites2 <- scaled_sites
#methyl_func_parallel(methyl = Y, sites = scaled_sites2, k_star = k_star,n_star =  n_star_list[[1]], X = X, ind_na_sub = NULL)

methyl_func_parallel(Y,scaled_sites,k_star,n_star_list[[1]], X,NULL)



methyl = Y; sites = scaled_sites; k_star = k_star;n_star =  n_star_list[[1]]; X = X; ind_na_sub = NULL
length(n_star_list[[1]])
dim(Y)


#---


# 1M - 90% N*
#             RMSE        R2 RMSE_sub R2_sub    time
# ols_gasp 0.1167959 0.8598460      NaN    NaN 109.427
# gasp     0.1074292 0.8814245      NaN    NaN  98.403
# null     0.1512797 0.7648680      NaN    NaN   1.129
#------------------------------------------------------------------------------
# all data; 99% N*
# RMSE        R2 RMSE_sub R2_sub    time
# ols_gasp 0.1462965 0.7349718      NaN    NaN 771.050
# gasp     0.1284976 0.7955373      NaN    NaN 732.916
# null     0.1388718 0.7611900      NaN    NaN  17.198
#-----------------------------------------------------------------------------
# all data; 90% N* ; K* all but 1 & 9; same as above too
# RMSE        R2 RMSE_sub R2_sub     time
# ols_gasp 0.1140599 0.8390043      NaN    NaN 1732.309
# gasp     0.1055445 0.8621460      NaN    NaN 1591.633
# null     0.1388924 0.7612711      NaN    NaN   17.396
#-----------------------------------------------------------------------------
# all data; 90% N* ; K* all but 6 rows chosen by selecting one row of possible combination
# RMSE        R2 RMSE_sub R2_sub     time
# ols_gasp 0.09973562 0.8765287      NaN    NaN 1840.185
# gasp     0.08447753 0.9114175      NaN    NaN 1599.213
# null     0.11277985 0.8421195      NaN    NaN   16.692
#-----------------------------------------------------------------------------
# same as previous but with no covariates (i passed a single dummy covariate)
#----------------------------------------------------------------------------------
# Male(24); Age only; alpha=0.05, min freq = 1, DMR only; chr7: K* (22) 1&3 are for training.
# Scenario 4
# RMSE        R2 RMSE_sub R2_sub    time
# ols_gasp 0.1024435 0.8559119      NaN    NaN 835.626
# gasp     0.1015363 0.8584526      NaN    NaN 718.261
# null     0.1084146 0.8386255      NaN    NaN   8.300
#--------------------------------------------------------------------------------------
# as the previous but N* is 90% random + DMR (scenario 3)
# RMSE        R2   RMSE_sub    R2_sub    time
# ols_gasp 0.09709191 0.8643102 0.10237422 0.8561066 310.715
# gasp     0.09389646 0.8730947 0.09874733 0.8661217 305.492
# null     0.10709643 0.8349060 0.10841457 0.8386255   8.296
#-------------------------------------------------------------------------------------------
# as previous but N^ is 90% random only (scenario 1)
# RMSE        R2   RMSE_sub    R2_sub    time
# ols_gasp 0.09058889 0.8815359 0.08654074 0.8971745 441.181
# gasp     0.08628881 0.8925155 0.08240335 0.9067714 386.222
# null     0.10701049 0.8346936 0.10285006 0.8547659   8.373
#----------------------------------------------------------------------------------
run_model_on_chr("chr12", subset=NA)  