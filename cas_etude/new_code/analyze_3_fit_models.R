run_model_on_chromosome <- function(chromosome, Age.Only=TRUE, Male.Only=TRUE,alpha=.05,
                             min_freq=1, min_k=2,subset=NA, no_cores=9,
                             floor_by=1e7, middle_point=FALSE, note = "", correction = function(alpha, N)(alpha/N)){
    
  setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")  
  Y = readRDS(paste0("new_data/Ydat_common_",chromosome, note, ".rds"))
  sites = readRDS(paste0("new_data/sites_common_",chromosome, note, ".rds"))
  X = readRDS(paste0("new_data/Xdat_common_",chromosome, note, ".rds")) %>%
    as.data.frame() %>% 
    mutate(AGE = (AGE - mean(AGE))/ sd(AGE) ) %>% 
    mutate(FEMALE = 1 - MALE) %>%
    mutate(MALE = (MALE - mean(MALE))/ sd(MALE)) %>%
    as.matrix()
  p_values = readRDS(paste0("new_data/p_values_",chromosome, note, ".rds"))

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
   
  corrected_alpha = correction(alpha, N)
  dmr_regions = get_dmr_regions(p_values, sites, corrected_alpha, floor_by = floor_by,
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
    X = as.matrix(X[,c("AGE", "MALE")]) 
    #cols_to_drop <- apply(X, 2, function(x) var(x, na.rm = TRUE) == 0)
    #X <- as.matrix(X[, !cols_to_drop])
    print(dim(X))
    print(X[1:3,])
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
                                                 Age.Only,"_",alpha,"x",min_freq, note, ".Rdata"))
  return(res) 
}
#-----------------------------------------------------------------------------------------------
# Example run
# res = run_model_on_chr("chr17", subset=NA, min_freq = 1, Age.Only = TRUE,
#                        no_cores = 5, middle_point=TRUE, floor_by=1e6,
#                        alpha=.2, min_k = 4)
