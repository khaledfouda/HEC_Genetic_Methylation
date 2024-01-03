methyl_func = function(methyl,sites,k_star,n_star,X, ind_na_sub = NULL){
  
  Y_obs = methyl
  Y_obs[k_star,n_star] = NA
  
  Y = methyl
  
  
  time_ols_gasp = system.time({
    obj_ols_gasp = test_ols_gasp(Y_obs, sites, X)
    obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)
    Y_ols_gasp = pred_fgasp(obj_ols_gasp)
    Y_ols_gasp[Y_ols_gasp < 0] = 0
    Y_ols_gasp[Y_ols_gasp > 1] = 1
  })[3]
  
  
  
  
  time_gasp = system.time({
    obj_gasp =  fit_fgasp(Y_obs, sites)
    Y_gasp = pred_fgasp(obj_gasp)
    Y_gasp[Y_gasp < 0] = 0
    Y_gasp[Y_gasp > 1] = 1
  })[3]
  
  
  
  time_null = system.time({Y_null = apply(Y_obs,2,function(x) {
    x[is.na(x)] = mean(x,na.rm=T) 
    return(x)
  })})[3]
  
  ind_na = is.na(Y_obs)
  Y_mean = apply(Y,2,mean)
  
  tab_res = data.frame(RMSE = c(RMSE(Y[ind_na], Y_ols_gasp[ind_na]),
                                RMSE(Y[ind_na], Y_gasp[ind_na]),
                                RMSE(Y[ind_na], Y_null[ind_na])),
                       R2 = c(R2(Y[ind_na], Y_ols_gasp[ind_na]),
                                R2(Y[ind_na], Y_gasp[ind_na]),
                                R2(Y[ind_na], Y_null[ind_na])),
                       RMSE_sub = c(RMSE(Y[k_star,ind_na_sub], Y_ols_gasp[k_star,ind_na_sub]),
                                   RMSE(Y[k_star,ind_na_sub], Y_gasp[k_star,ind_na_sub]),
                                   RMSE(Y[k_star,ind_na_sub], Y_null[k_star,ind_na_sub])),
                       R2_sub = c(R2(Y[k_star,ind_na_sub], Y_ols_gasp[k_star,ind_na_sub]),
                                    R2(Y[k_star,ind_na_sub], Y_gasp[k_star,ind_na_sub]),
                                  R2(Y[k_star,ind_na_sub], Y_null[k_star,ind_na_sub])),
                       time = c(time_ols_gasp,
                                time_gasp,
                                time_null)
                       )
  rownames(tab_res) = c("ols_gasp","gasp","null")
  
  return(tab_res)
}