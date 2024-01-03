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

load(file = "data/res_mcci4.Rdata")

for (seed in 51:n_rep){
  print(seed)
  K = N.popu = 600 # Population size, the value of $N$ in the paper.
  N = N.ques = 1200  # Number of questions, the value of $L$ in the paper.
  D = N.covariates = 20
  n.low.rank = 10 # Rank of the response matrix $B_N^*$ in the paper.
  set.seed(seed)
  beta = matrix(rnorm(N.ques * N.covariates, mean = 0, sd = 1), nrow = N.covariates)  # Generate the regression parameter matrix $\beta_N$
  beta[beta < quantile(beta, 0.5)] = 0  # Smaller half values of $\b_N^*$ are set to 0.
  set.seed(seed)
  X = X.cov = matrix(rnorm(N.popu * N.covariates, mean = 0, sd = 1), nrow = N.popu)  # Generate the covariate matrix $X_N$.
  set.seed(seed)
  BL = matrix(rnorm(N.popu * n.low.rank, mean = 0, sd = 1), nrow = N.popu)  # Generate  $B_L$
  set.seed(seed)
  BR = matrix(rnorm(N.ques * n.low.rank, mean = 0, sd = 1), nrow = N.ques)  # Generate $B_R$
  PX = X.cov %*% solve(t(X.cov) %*% X.cov) %*% t(X.cov)  # Generate the projection matrix to guarantee the orthogonality of $X_N$ and $B_N^*$ on the population level.
  Eye = diag(1, dim(X.cov)[1])  # Generate the identity matrix
  PXp = Eye - PX
  B0 = PXp %*% BL %*% t(BR)
  Y = X.cov %*% beta + B0
  sig=sqrt(var(as.vector(Y)))
  set.seed(seed)
  eps = matrix(rnorm(K*N,0,sig),K,N)
  
  
  
  
  Y_obs = Y + eps
  
  
  # gamma.i.1 = matrix(rnorm(3 * K, -0.1, 0.1), nrow = 3)  # Generate $\gamma_j$ for $j=1,\ldots,3$.
  # gamma.i.0 = matrix(rnorm(1 * K, 0.3, 0.1), nrow = 1)  # Generate $\gamma_0$
  # mp = 1/(1 + exp(-(cbind(1, X.cov[, 1:3]) %*% rbind(gamma.i.0, gamma.i.1))))  # Only the first three components are used for the response model.
  # omega=matrix(rbinom(K*N,1,mp),K,N)
  # Y_obs[omega==0] = NA
  
  temp = Y_obs
  set.seed(seed)
  random_NA = sample(1:length(Y),round(length(Y)*0.9))
  temp[random_NA] = NA
  col_with_value = unique(round(seq(1,ncol(Y),length.out = round(N/2))))
  temp[,col_with_value] = Y_obs[,col_with_value]
  Y_obs = temp
  sum(is.na(Y_obs))/(K*N)
  
  apply(is.na(Y_obs),1,sum)
  
  # n_star = sort(sample(1:N,round(0.9*N)))
  # k_star = 11:K
  # Y_obs[k_star,n_star] = NA
  # sum(is.na(Y_obs))/(K*N)
  
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
      while(sum(is.na(Y_ols_gasp))>0  & count < 10){
        print(count)
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
  
  
  Y_null = apply(Y_obs,2,function(x) {
    x[is.na(x)] = mean(x,na.rm=T) 
    return(x)
  })
  
  # ind_na = is.na(Y_obs)
  # cat("ols_gaSP",MAE(Y[ind_na],Y_ols_gasp[ind_na]),
  #     " mcci",MAE(Y[ind_na], Y_mcci[ind_na]),
  #     " null",MAE(Y[ind_na], Y_null[ind_na]),
  #     "\n")
  
  
  norm(Y_ols_gasp-Y,type="F")/sqrt(K*N)
  norm(Y_mcci-Y,type="F")/sqrt(K*N)
  norm(Y_null-Y,type="F")/sqrt(K*N)
  
  
  ind_na = is.na(Y_obs)
  cat("ols_gaSP",RMSE(Y[ind_na],Y_ols_gasp[ind_na])," mcci",RMSE(Y[ind_na], Y_mcci[ind_na]),"\n")
  
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
    save(res,file = "data/res_mcci4.Rdata")
  }
  
}

save(res,file = "data/res_mcci4.Rdata")