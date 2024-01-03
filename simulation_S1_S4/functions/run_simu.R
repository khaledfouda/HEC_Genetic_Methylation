run_simu = function(N,
                    K,
                    D,
                    n_star,
                    k_star,
                    method = c(1, 2),
                    sX = rep(1, length(sites)),
                    seed = as.numeric(Sys.time()),
                    discrete = FALSE,
                    mispe = F,
                    wrong_X = F,
                    sinuzo = F,
                    V = NULL,
                    trend = F) {
  
  
  set.seed(seed)
  sites = sort(runif(N))
  set.seed(seed)
  beta = c(runif(D,2,5),runif(K-D,10,1000))
  set.seed(seed)
  nugget = c(rep(0,D),runif(K-D,0.001,0.05))

  simu = simu_GaSP(D,sites,beta, nugget,0,seed,discrete,mispe,sinuzo,V_save = V,trend = trend)
  
  
  Y = simu$Y
  
  
  
  X = simu$X
  
  
  #X
  
  # a = which.min(X)
  # b = which.max(X)
  # 
  # plot(sites,Y[a,], cex = 0.5, pch = 16, col="blue")
  # lines(sites,X[a,]%*%simu$B,col="blue")
  # points(sites,Y[b,], cex = 0.5, pch = 16, col="red")
  # lines(sites,X[b,]%*%simu$B,col="red")
  
  #X = runif(K,-1,1)
  
  if ( wrong_X){
    if (discrete){
      if(D!=2){
        stop("in discrete case only D=2 allowed")
      }
      x = sample(rep(c(1,0),ceiling(K/2))[1:K])
      X = matrix(0,ncol = 2, nrow= K)
      X[x==1,1] = 1
      X[x==0,2] = 1
    } else {
      X = matrix(runif(K*D),ncol=D)
    }
  }
  
  
  
  
  par(mfrow=c(1,1))
  # par(mar=c(4.1,2.1,4.1,2.1))
  # plot(sites,simu$B[1,],type="l",ylim=range(c(-simu$B,simu$B)), ylab = "Xb(s)", main = "Moyenne",col=3)
  # lines(sites,simu$B[2,],type="l")


  # plot(sites,Y[1,], ylim = range(Y),pch = 16,type = "l", cex=0.5,col=X[1]+2, ylab = "Y", main = "Simulation")
  # points(sites[-n_star],Y[1,-n_star], ylim = range(Y),pch = 16, cex = 0.5)
  # for ( i in 2:10){
  #   points(sites,Y[i,], ylim = range(Y),type = "l",pch = 16, cex = 0.5,col=X[i]+2)
  #   points(sites[-n_star],Y[i,-n_star], ylim = range(Y),pch = 16, cex = 0.5)
  # }
  # 
  # 
  # legend("bottomright",c("X=1","X=-1"),lty=1,col=c(3,1))
  # 

  ## missing values
  Y_obs = Y

  Y_obs[k_star,n_star] = NA
  
  
  
  
  #X = cbind(X,rbinom(K,1,0.5),runif(K),rnorm(K))
  
  if(1 %in% method){
    t_ols_gasp = system.time({
      count = 0
      test = try({
        obj_ols_gasp = test_ols_gasp(Y_obs,sites,X, tol = 1e-6)
        obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)#,B = simu$B[-n_star])#, A = simu$A)#,A_output_t = rbind(simu$B,simu$V)[,-n_star])
        Y_ols_gasp = pred_ols_gasp(obj_ols_gasp)
        while(sum(is.na(Y_ols_gasp))>0 & count < 20){
          obj_ols_gasp = test_ols_gasp(Y_obs,sites,X)
          obj_ols_gasp = fit_ols_gasp(obj_ols_gasp)#,B = simu$B[-n_star])#, A = simu$A)#,A_output_t = rbind(simu$B,simu$V)[,-n_star])
          Y_ols_gasp = pred_ols_gasp(obj_ols_gasp)
          count = count + 1
        }
      })
    })
  }
  
  if(class(test)[1] == "try-error" | !1 %in% method){
    Y_ols_gasp = NULL
    t_ols_gasp = rep(NA,3)
  }
  
  if(2 %in% method){
    t_gasp  = system.time({
      test = try({
        obj_gasp =  fit_fgasp(Y_obs,sites)#,A = cbind(X,simu$A), A_output_t = rbind(simu$B,simu$V)[,-n_star])
        Y_gasp = pred_fgasp(obj_gasp)
        while(sum(is.na(Y_gasp))>0 & count < 20){
          obj_gasp =  fit_fgasp(Y_obs,sites)#,A = cbind(X,simu$A), A_output_t = rbind(simu$B,simu$V)[,-n_star])
          Y_gasp = pred_fgasp(obj_gasp)
          count = count + 1
        }
      })})
  }
  if(class(test)[1] == "try-error" | !2 %in% method){
    Y_gasp = NULL
    t_gasp = rep(NA,3)
  }
  
  ind_na = is.na(Y_obs)
  cat("ols_gaSP",RMSE(Y[ind_na],Y_ols_gasp[ind_na])," gaSP",RMSE(Y[ind_na],Y_gasp[ind_na]),"\n")
  
  
  ### essai du code du mec
  # Aobs = Y_obs
  # Aobs[ind_na] = 0
  # diagD = rep(1,nrow(Y_obs))
  # sample.mis.SMC = SMCfit(Aobs = Aobs, X = X, diagD, 0.5,0.5, 1)
  # fitted.beta.sample.SMC = sample.mis.SMC$beta
  # fitted.B.sample.SMC = sqrt(diagD) * sample.mis.SMC$Bhat
  # fitted.sample.SMC = X %*% sample.mis.SMC$beta + sqrt(diagD) * sample.mis.SMC$Bhat
  # 
  # RMSE(Y[ind_na],fitted.sample.SMC[ind_na])
  # 
  # 
  
  if(3 %in% method){
    t_gls_gasp = system.time({
      test = try({
        Y_gls_gasp = gls_gasp(formula, Y_obs, sites, X, D = D)
      })
    })
  }

  if(class(test)[1] == "try-error" | !3 %in% method){
    Y_gls_gasp = NULL
    t_gls_gasp = rep(NA,3)
  }


  if(4 %in% method){
    t_gls_gasp_up = system.time({
      test = try({
        Y_gls_gasp_up = gls_up_gasp(formula, Y_obs, sites, X, D = D)
      })})
  }

  if(class(test)[1] == "try-error" | !4 %in% method){
    Y_gls_gasp_up = NULL
    t_gls_gasp_up = rep(NA,3)
  }
  
  
  return(
    list(
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
      Y_gasp = Y_gasp,
      Y_gls_gasp = Y_gls_gasp,
      Y_gls_gasp_up = Y_gls_gasp_up,
      time = c(t_ols_gasp[3], t_gasp[3], t_gls_gasp[3], t_gls_gasp_up[3])
    )
  )
  

}

