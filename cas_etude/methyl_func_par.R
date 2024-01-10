methyl_func_parallel <- function(methyl,sites,k_star,n_star,X, ind_na_sub = NULL) {

  
  
  #----------------------------------------------------------
  Y_obs = methyl
  Y_obs[k_star, n_star] = NA
  Y = methyl
  ind_na = is.na(Y_obs)
  Y_mean = apply(Y, 2, mean)
  #-----------------------------------------------------------------------
  tasks = list(
    ols_gasp = function() {
      time = system.time({
        obj_ols_gasp = test_ols_gasp(Y_obs, sites, X)
        obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)
        Y_ols_gasp = pred_fgasp(obj_ols_gasp)
        Y_ols_gasp[Y_ols_gasp < 0] = 0
        Y_ols_gasp[Y_ols_gasp > 1] = 1
      })[3]
      list(Y = Y_ols_gasp, time = time)
    },
    gasp = function() {
      time = system.time({
        obj_gasp =  fit_fgasp(Y_obs, sites)
        Y_gasp = pred_fgasp(obj_gasp)
        Y_gasp[Y_gasp < 0] = 0
        Y_gasp[Y_gasp > 1] = 1
      })[3]
      list(Y = Y_gasp, time = time)
    },
    null = function() {
      time = system.time({Y_null = apply(Y_obs,2,function(x) {
        x[is.na(x)] = mean(x,na.rm=T) 
        return(x)
      })})[3]
      list(Y = Y_null, time =time)
    }
  )
  #-----------------------------------------------------
  inner_cl <- makeCluster(3)
  registerDoParallel(inner_cl)
  nothing = clusterEvalQ(inner_cl, {
    require(Jmisc)
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
    #source("methyl_func.R")
    sourceAll(path="../functions/")
  })
  #print("hi")
  clusterExport(inner_cl, c("Y_obs", "Y", "sites", "X", "ind_na_sub", "tasks", "ind_na", "Y_mean","scaled_sites"))
                            #"test_ols_gasp","fit_ols_gasp", "pred_fgasp", "fit_fgasp", 
  #print("hi2")
  results <- foreach(task = iter(tasks)) %dopar% {
    task()
  }
  #results = parLapply(inner_cl, tasks, function(task) task())  
  #print("hi3")
  stopCluster(inner_cl)
  #------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------
  tab_res = data.frame(
    RMSE = c(RMSE(Y[ind_na], results[[1]]$Y[ind_na]),
             RMSE(Y[ind_na], results[[2]]$Y[ind_na]),
             RMSE(Y[ind_na], results[[3]]$Y[ind_na])),
    R2 = c(R2(Y[ind_na], results[[1]]$Y[ind_na]),
           R2(Y[ind_na], results[[2]]$Y[ind_na]),
           R2(Y[ind_na], results[[3]]$Y[ind_na])),
    RMSE_sub = c(RMSE(Y[k_star, ind_na_sub], results[[1]]$Y[k_star, ind_na_sub]),
                 RMSE(Y[k_star, ind_na_sub], results[[2]]$Y[k_star, ind_na_sub]),
                 RMSE(Y[k_star, ind_na_sub], results[[3]]$Y[k_star, ind_na_sub])),
    R2_sub = c(R2(Y[k_star, ind_na_sub], results[[1]]$Y[k_star, ind_na_sub]),
               R2(Y[k_star, ind_na_sub], results[[2]]$Y[k_star, ind_na_sub]),
               R2(Y[k_star, ind_na_sub], results[[3]]$Y[k_star, ind_na_sub])),
    time = c(results[[1]]$time, results[[2]]$time, results[[3]]$time)
  )
  
  rownames(tab_res) = c("ols_gasp", "gasp", "null")
  
  return(tab_res)
}


system.time({results=methyl_func_parallel(Y,scaled_sites,k_star,n_star_list[[1]], X,NULL)})
