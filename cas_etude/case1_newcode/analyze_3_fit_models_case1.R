run_model_case1 <- function(alpha = .05,
                            no_cores = 9,
                            floor_by = 1e3,
                            note = "",
                            correction = FALSE) {
  setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")
  
  
  source("methyl_func.R")
  sourceAll(path = "../functions/")
  
  ## get data from fastgasp article
  load("data/Methylation_level_dense_Chr1.Rda")
  Y = t(seq.data.chr1)
  selec = c(1:2, 4:9, 12:14, 16:19)
  Y = Y[selec,]
  colnames(Y) = rownames(Y) = NULL
  X1 = c(rep("Gifford", 8), rep("roadmap1", 7))
  #X = c(rep(1, 8), rep(0, 7))
  #X = matrix((X - mean(X)) / sd(X), ncol=1) 
  X = fact2mat(X1)
  sites = Index_seq
  p_values = readRDS(paste0("new_data/case1_p_values", ".rds"))
  N = ncol(Y)
  K = nrow(Y)
  print("Starting..")
  
  print(paste0("Pre-correction alpha = ",alpha))
  
  if(correction) alpha = alpha / N
  print(paste0("Post-correction alpha = ",alpha))
  
  dmr_regions <- get_dmr_regions_case1(p_values,
                                       sites,
                                       alpha,
                                       floor_by)
  
  print("Starting..")
  
  
  ind_dmr = sort(which(sites %in% dmr_regions))
  sites = scale_01(sites)
  
  
  
  k_star = (1:K)[-c(1, 9)]
  
  #---------------
  n_star_list = list(sort(unique(c(round(
    seq(1, N, length.out = round(N * 0.9))
  )))), #scenario 1,
  sort(unique(c(
    round(seq(1, N, length.out = round(N * 0.9))), ind_dmr
  ))), #scenario 3
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
    sourceAll(path = "../functions/")
  })
  print(dim(Y))
  print(length(n_star_list[[1]]))
  print(length(n_star_list[[2]]))
  print(length(n_star_list[[3]]))
  print(dim(X))
  print(length(k_star))
  
  
  results <- foreach(i = 1:9, .combine = "rbind") %dopar% {
    if (i %in% 1:3) {
      ind_na_sub = NULL
      n_star <- n_star_list[[1]]
    } else if (i %in% 4:6) {
      ind_na_sub = ind_dmr
      n_star <- n_star_list[[2]]
    } else if (i %in% 7:9) {
      ind_na_sub = ind_dmr
      n_star <- n_star_list[[3]]
    }
    
    methyl <- Y
    methyl[k_star, n_star] <- NA
    
    if (i %in% c(1, 4, 7)) {
      ind_na = is.na(methyl)
      time = system.time({
        obj_ols_gasp = test_ols_gasp(methyl, sites, X)
        obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)
        Y_pred = pred_fgasp(obj_ols_gasp)
        Y_pred[Y_pred < 0] = 0
        Y_pred[Y_pred > 1] = 1
      })[3]
      model = "LMCC"
      
    } else if (i %in% c(2, 5, 8)) {
      ind_na = is.na(methyl)
      time = system.time({
        obj_gasp =  fit_fgasp(methyl, sites)
        Y_pred = pred_fgasp(obj_gasp)
        Y_pred[Y_pred < 0] = 0
        Y_pred[Y_pred > 1] = 1
      })[3]
      model = "GASP"
      
    } else if (i %in% c(3, 6, 9)) {
      ind_na = is.na(methyl)
      time = system.time({
        Y_pred = apply(methyl, 2, function(x) {
          x[is.na(x)] = mean(x, na.rm = T)
          return(x)
        })
      })[3]
      model = "Naive"
    }
    out = data.frame(
      RMSE = RMSE(Y[ind_na],  Y_pred[ind_na]),
      R2 =  R2(Y[ind_na], Y_pred[ind_na]),
      RMSE_dmr = RMSE(Y[k_star, ind_na_sub],  Y_pred[k_star, ind_na_sub]),
      R2_dmr =  R2(Y[k_star, ind_na_sub],  Y_pred[k_star, ind_na_sub]),
      time = time,
      model = model,
      n_star = length(n_star),
      miss_all = round(sum(ind_na) / (N * K), 3)
    )
    
    if (i %in% 1:3) {
      out$scenario = "1"
    } else if (i %in% 4:6) {
      out$scenario = "3"
    } else if (i %in% 7:9) {
      out$scenario = "4"
    }
    
    out
  }
  
  
  stopCluster(cl)
  
  res = as.data.frame(results) %>%
    arrange(scenario) %>%
    mutate(k_star = length(k_star),
           N = N,
           K = K) %>%
    mutate(miss_r = round(n_star / N, 2),
           miss_c = round(k_star / K, 2))
  row.names(res) <- NULL
  print(res)
  save(
    res,
    file = paste0(
      "results/res_case1_",
      alpha,
      note,
      ".Rdata"
    )
  )
  return(res)
}
#-----------------------------------------------------------------------------------------------
# Example run
# res = run_model_on_chr("chr17", subset=NA, min_freq = 1, Age.Only = TRUE,
#                        no_cores = 5, middle_point=TRUE, floor_by=1e6,
#                        alpha=.2, min_k = 4)
