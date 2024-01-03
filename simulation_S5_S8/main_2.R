source("MCCI_functions.R")

library(mvmeta)
library(Jmisc)
library(Rcpp)


## functions to estimate GP
sourceAll("../functions/")

sourceAll("../simulations/functions/")

n_rep = 100

res = list()
length(res) = n_rep



for (seed in 1:n_rep){
  print(seed)
  K = N.popu = 100  # Population size, the value of $N$ in the paper.
  N = N.ques = 1000  # Number of questions, the value of $L$ in the paper.
  D = N.covariates = 20
  n.low.rank = 10 # Rank of the response matrix $B_N^*$ in the paper.
  set.seed(seed)
  beta = matrix(rnorm(N.ques * N.covariates, mean = 0, sd = 1), nrow = N.covariates)  # Generate the regression parameter matrix $\beta_N$
  beta[beta < quantile(beta, 0.5)] = 0  # Smaller half values of $\b_N^*$ are set to 0.
  X = X.cov = matrix(rnorm(N.popu * N.covariates, mean = 0, sd = 1), nrow = N.popu)  # Generate the covariate matrix $X_N$.
  BL = matrix(rnorm(N.popu * n.low.rank, mean = 0, sd = 1), nrow = N.popu)  # Generate  $B_L$
  BR = matrix(rnorm(N.ques * n.low.rank, mean = 0, sd = 1), nrow = N.ques)  # Generate $B_R$
  PX = X.cov %*% solve(t(X.cov) %*% X.cov) %*% t(X.cov)  # Generate the projection matrix to guarantee the orthogonality of $X_N$ and $B_N^*$ on the population level.
  Eye = diag(1, dim(X.cov)[1])  # Generate the identity matrix
  PXp = Eye - PX
  B0 = PXp %*% BL %*% t(BR)
  Y = X.cov %*% beta + B0
  sig=sqrt(var(as.vector(Y)))
  eps = matrix(rnorm(K*N,0,sig),K,N)
  
  
  Y_obs = Y + eps
  
  n_star = sort(sample(1:N,round(0.9*N)))
  k_star = 11:K
  Y_obs[k_star,n_star] = NA
  
  
  
  #Y_obs[sample(1:length(Y_obs),N*0.4)] = NA
  
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
        mcci_fit = MCCIfit(Aobs,X,nfolds=5,SVT_type="FRSVT",lambda1_grid=seq(0,1000,length=30),
                           lambda2_gridlength=50,alpha_grid=1)
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
    sites = seq(0,1,length.out=ncol(Y_obs))
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
  cat("ols_gaSP",MAE(Y[ind_na],Y_ols_gasp[ind_na])," mcci",MAE(Y[ind_na], Y_mcci[ind_na]),"\n")
  
  res[[seed]] = list(
    X = X,
    A = BL,
    V = BR,
    beta = beta,
    Y = Y,
    sites = sites,
    Y_obs = Y_obs,
    Y_ols_gasp = Y_ols_gasp,
    Y_mcci = Y_mcci,
    time = c(t_ols_gasp[3], t_mcci[3])
  )
  
  if(seed %%10 == 0){
    save(res,file = "data/res_mcci2.Rdata")
  }
  
}

save(res,file = "data/res_mcci2.Rdata")