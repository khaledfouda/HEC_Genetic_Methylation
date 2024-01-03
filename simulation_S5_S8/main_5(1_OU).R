source("MCCI_functions.R")

library(mvmeta)
library(Jmisc)
library(Rcpp)


## functions to estimate GP
sourceAll("../functions/")

sourceAll("../simulation_S1_S4_OU/functions/")

n_rep = 100

res = list()
length(res) = n_rep

for (seed in 1:n_rep){
  print(seed)
  ## parameters ##
  N = 1000
  K = 100
  D = 20
  k_star = c(11:K)
  n_star = sort(sample(1:N,round(N*0.9)))
  
  ## GP simulations ##
  set.seed(seed)
  sites = sort(runif(N))
  set.seed(seed)
  beta = c(runif(D,2,5),runif(K-D,10,1000))
  set.seed(seed)
  nugget = c(rep(0,D),runif(K-D,0.001,0.05))
  simu = simu_GaSP(D,sites,beta, nugget,0,seed)
  
  
  Y = simu$Y
  X = simu$X
  
  
  Y_obs = Y
  
  Y_obs[k_star,n_star] = NA
  sum(is.na(Y_obs))/(K*N)
  
  Aobs = Y_obs
  Aobs[is.na(Aobs)] = 0
  
  t_mcci = system.time({
    count = 0
    test = try({
      mcci_fit = MCCIfit(Aobs,X, lambda2_grid = 50)
      Y_mcci = mcci_fit$Ahat
    })
    while(class(test)[1] == "try-error" & count < 2){
      test = try({
        mcci_fit = MCCIfit(Aobs,X, lambda2_grid = round(50/(count+2)))
        Y_mcci = mcci_fit$Ahat
      })
      count = count + 1
      print(count)
    }
  })
  
  if(class(test)[1] == "try-error" ){
    Y_mcci = NULL
    t_mcci = rep(NA,3)
  }
  
  
  t_ols_gasp = system.time({
    count = 0
    sites = sort(runif(N))
    test = try({
      obj_ols_gasp = test_ols_gasp(Y_obs,sites,X)
      obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)#,B = simu$B[-n_star])#, A = simu$A)#,A_output_t = rbind(simu$B,simu$V)[,-n_star])
      Y_ols_gasp = pred_ols_gasp(obj_ols_gasp)
      while(sum(is.na(Y_ols_gasp))>0 & count < 10){
        obj_ols_gasp = test_ols_gasp(Y_obs,sites,X)
        obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)#,B = simu$B[-n_star])#, A = simu$A)#,A_output_t = rbind(simu$B,simu$V)[,-n_star])
        Y_ols_gasp = pred_ols_gasp(obj_ols_gasp)
        count = count + 1
      }
    })
  })
  
  if(class(test)[1] == "try-error"){
    Y_ols_gasp = NULL
    t_ols_gasp = rep(NA,3)
  }
  
  
  
  ind_na = is.na(Y_obs)
  cat("ols_gaSP",RMSE(Y[ind_na],Y_ols_gasp[ind_na])," mcci",RMSE(Y[ind_na], Y_mcci[ind_na]),"\n")
  #
  
  res[[seed]] = list(
    X = X,
    A = simu$A,
    B = simu$B,
    V = simu$V,
    beta = beta,
    nugget = nugget,
    Y = Y,
    sites = sites,
    Y_obs = Y_obs,
    Y_ols_gasp = Y_ols_gasp,
    Y_mcci = Y_mcci,
    time = c(t_ols_gasp[3], t_mcci[3])
  )
  if(seed %%10 == 0){
    save(res,file = "data/res_mcci5.Rdata")
  }
  
}

save(res,file = "data/res_mcci5.Rdata")