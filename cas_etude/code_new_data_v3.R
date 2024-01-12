


get_dmr_regions <- function(p_values, sites, N, alpha=0.1, floor_by=1e7, min_freq=2, return_seq=FALSE,
                            middle_point=FALSE){
  if(middle_point == TRUE){
    loci = sort(sites[which(p_values < (alpha/N))])
    if(return_seq == TRUE){
      sequen = c()
      step_size = round(floor_by/2)
      end_site = 1
      for(i in 1:length(loci)){
        start_site = max( end_site, loci[i]-step_size)
        end_site = loci[i] + step_size
        sequen = c(sequen, start_site: end_site)
      }
      return(sequen)
    }
    return(loci)
  }
  
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



#chromosome="chr12"; Age.Only=TRUE; Male.Only=TRUE;alpha=.05; min_freq=1; min_k=2;subset=1e5
#-----------------------
run_model_on_chr <- function(chromosome, Age.Only=TRUE, Male.Only=TRUE,alpha=.05,
                             min_freq=1, min_k=2,subset=NA, no_cores=9,
                             floor_by=1e7, middle_point=FALSE){
    
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
  dmr_regions = get_dmr_regions(p_values, sites, N, alpha, floor_by = floor_by,
                                min_freq = min_freq, middle_point = middle_point, return_seq = T)
  ind_dmr = sort(which(sites %in% dmr_regions))
  sites = scale_01(sites)
  
  
  k_star <-
    (X %>%
    as.data.frame() %>%
    mutate(index = 1:n()) %>%  
    group_by(MALE, AML, APL, BONE_MARROW) %>% 
    filter(row_number() != 1) %>% 
    ungroup())$index  
  
  if( (K-length(k_star)) < min_k ){
    extra_k = sample(1:length(k_star), (min_k-K+length(k_star)),replace = FALSE )
    k_star <- k_star[-extra_k]
  }
  if(Age.Only == TRUE){
    X = as.matrix(X[,"AGE"])
  }else{
    # check for columns with single values and drop them
    cols_to_drop <- apply(X, 2, function(x) var(x, na.rm = TRUE) == 0)
    X <- as.matrix(X[, !cols_to_drop])
    print(dim(X))
  }
  
  #---------------
  n_star_list = list(sort(unique(c(round(seq(1,N,length.out=round(N*0.9)))))), #scenario 1,
                     sort(unique(c(round(seq(1,N,length.out=round(N*0.9))),ind_dmr))), #scenario 3
                     sort(ind_dmr)) #scenario 4
  
  
  
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  nothing = clusterEvalQ(cl, {
    require(Jmisc)
    require(DiceEval)
    require(mgcv)
    require(Rcpp)
    require(RColorBrewer)
    require(ggpubr)
    require(ZIprop)
    require(corrplot)
    require(tidyverse)
    require(magrittr)
    sourceAll(path="../functions/")
  })
  print(dim(Y))
  print(length(n_star_list[[1]]))
  print(length(n_star_list[[2]]))
  print(length(n_star_list[[3]]))
  print(dim(X))
  print(length(k_star))
  
  
  results <- foreach(i = 1:9, .combine = "rbind") %dopar% {
    if(i %in% 1:3){
      ind_na_sub = NULL
      n_star <- n_star_list[[1]]
    }else if(i %in% 4:6){
      ind_na_sub = ind_dmr
      n_star <- n_star_list[[2]]
    }else if(i %in% 7:9){
      ind_na_sub = ind_dmr
      n_star <- n_star_list[[3]]
    }
  
    methyl <- Y  
    methyl[k_star, n_star] <- NA
    
    if(i %in% c(1,4,7)){
      
      ind_na = is.na(methyl)
      time = system.time({
          obj_ols_gasp = test_ols_gasp(methyl, sites, X)
          obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)
          Y_pred = pred_fgasp(obj_ols_gasp)
          Y_pred[Y_pred < 0] = 0
          Y_pred[Y_pred > 1] = 1
        })[3]
        model = "LMCC"
        
      }else if(i %in% c(2,5,8)){
        ind_na = is.na(methyl)
        time = system.time({
          obj_gasp =  fit_fgasp(methyl, sites)
          Y_pred = pred_fgasp(obj_gasp)
          Y_pred[Y_pred < 0] = 0
          Y_pred[Y_pred > 1] = 1
        })[3]
        model = "GASP"

      }else if(i %in% c(3,6,9)){
        ind_na = is.na(methyl)
        time = system.time({Y_pred = apply(methyl,2,function(x) {
          x[is.na(x)] = mean(x,na.rm=T) 
          return(x)
        })})[3]
        model = "Naive"
      }
    out = data.frame(RMSE = RMSE(Y[ind_na],  Y_pred[ind_na]),R2 =  R2(Y[ind_na],Y_pred[ind_na]),
                     RMSE_dmr =RMSE(Y[k_star, ind_na_sub],  Y_pred[k_star, ind_na_sub]),
                     R2_dmr =  R2(Y[k_star, ind_na_sub],  Y_pred[k_star, ind_na_sub]),
                     time = time, model = model, n_star = length(n_star), miss_all = round(sum(ind_na)/(N*K),3))
    
    if(i %in% 1:3){
      out$scenario = "1"
    }else if(i %in% 4:6){
      out$scenario = "3"
    }else if(i %in% 7:9){
      out$scenario = "4"
    }
    
    out
  }
   

  stopCluster(cl)

  res = as.data.frame(results) %>%
    arrange(scenario) %>%  
           mutate(k_star = length(k_star),
           N = N, K = K) %>% 
    mutate(miss_r = round(n_star/N,2),
           miss_c = round(k_star/K,2))
  row.names(res) <- NULL
  print(res)
  save(res, file = paste0("results/res_dmr1_",chromosome,"_MaleOnly_",Male.Only,"_AgeOnly_",
                                                 Age.Only,"_",alpha,"x",min_freq,".Rdata"))
  return(res) 
}
#-----------------------------------------------------------------------------------------------

for(chr in paste0("chr",c(7,8,11,12)))
res = run_model_on_chr(chr, subset=NA, min_freq = 1, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                       alpha=.1)  

for(chr in paste0("chr",c(7,8,11,12,17)))
  res = run_model_on_chr(chr, subset=NA, min_freq = 1, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                         alpha=.2)



for(chr in paste0("chr",c(7,8,11,12,17)))
  res = run_model_on_chr(chr, subset=NA, min_freq = 1, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                         alpha=.05)


for(chr in paste0("chr",c(7,8,11,12)))
  res = run_model_on_chr(chr, subset=NA, Age.Only = FALSE, min_freq = 1, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                         alpha=.1)  

for(chr in paste0("chr",c(7,8,11,12,17)))
  res = run_model_on_chr(chr, subset=NA, Age.Only = FALSE, min_freq = 1, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                         alpha=.2)



for(chr in paste0("chr",c(7,8,11,12,17)))
  res = run_model_on_chr(chr, subset=NA, min_freq = 1, Age.Only = FALSE, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                         alpha=.05)



res = run_model_on_chr("chr17", subset=NA, min_freq = 1, Age.Only = TRUE,
                       no_cores = 5, middle_point=TRUE, floor_by=1e6,
                       alpha=.2, min_k = 4)


#---------------------------------------------------------------------------------------------------
# alpha = .05
# min_frq = 1
# load(paste0("results/res_dmr_",chromosome,"_MaleOnly_TRUE_AgeOnly_TRUE_",
#                              alpha,"x",min_freq,".Rdata"))
# 
# 
# 
# res %>%
#   mutate(RMSE = round(RMSE, 3), R2 = round(R2, 3), RMSE_dmr = round(RMSE_dmr,3),
#          R2_dmr = round(R2_dmr, 3), time_min = round(time/60,1), time=NULL,
#          N = paste0(round(N/1e6,1),"M"),n_star = paste0(round(n_star/1e6,1),"M"))

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
# RMSE        R2  RMSE_dmr    R2_dmr    time model n_star miss_all scenario k_star      N  K miss_r miss_c
# 1 0.09760188 0.8699971       NaN       NaN 396.336  LMCC 674712    0.825        1     22 749680 24   0.90   0.92
# 2 0.09327078 0.8812789       NaN       NaN 392.516  GASP 674712    0.825        1     22 749680 24   0.90   0.92
# 3 0.11667007 0.8142386       NaN       NaN   7.424 Naive 674712    0.825        1     22 749680 24   0.90   0.92
# 4 0.10164926 0.8594600 0.1032892 0.8607157 217.862  LMCC 730866    0.894        3     22 749680 24   0.97   0.92
# 5 0.09871203 0.8674647 0.1006033 0.8678655 211.302  GASP 730866    0.894        3     22 749680 24   0.97   0.92
# 6 0.11673431 0.8146517 0.1170906 0.8210069   5.916 Naive 730866    0.894        3     22 749680 24   0.97   0.92
# 7 0.10329761 0.8606931 0.1032976 0.8606931 507.933  LMCC 561531    0.687        4     22 749680 24   0.75   0.92
# 8 0.10254534 0.8627148 0.1025453 0.8627148 434.452  GASP 561531    0.687        4     22 749680 24   0.75   0.92
# 9 0.11709064 0.8210069 0.1170906 0.8210069   5.699 Naive 561531    0.687        4     22 749680 24   0.75   0.92
#--------------

# chr 17; .01; middle; 1e6;
# RMSE        R2  RMSE_dmr    R2_dmr    time model n_star miss_all scenario k_star      N
# 1 0.09626129 0.9007937       NaN       NaN 346.638  LMCC 660405    0.825        1     22 733783
# 2 0.09541319 0.9025341       NaN       NaN 344.696  GASP 660405    0.825        1     22 733783
# 3 0.11958387 0.8468980       NaN       NaN   7.652 Naive 660405    0.825        1     22 733783
# 4 0.10647580 0.8793266 0.1113152 0.8790600 242.696  LMCC 703009    0.878        3     22 733783
# 5 0.10451162 0.8837377 0.1105054 0.8808132 235.554  GASP 703009    0.878        3     22 733783
# 6 0.11969288 0.8475083 0.1209823 0.8571418   5.731 Naive 703009    0.878        3     22 733783
# 7 0.11135368 0.8789763 0.1113537 0.8789763 824.619  LMCC 426050    0.532        4     22 733783
# 8 0.11358781 0.8740713 0.1135878 0.8740713 604.060  GASP 426050    0.532        4     22 733783
# 9 0.12098234 0.8571418 0.1209823 0.8571418   5.613 Naive 426050    0.532        4     22 733783

 