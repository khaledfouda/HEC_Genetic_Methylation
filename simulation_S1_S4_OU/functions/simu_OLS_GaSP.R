simu_OLS_GaSP= function(sites,A,beta, nugget,X,B, sX = rep(1,length(sites)), seed = as.numeric(Sys.time())){
  
  N = length(sites)
  K = nrow(A)
  D = nrow(A)

  
  ## random effect simulation
  set.seed(seed)
  simu = simu_GaSP(A, sites, beta, nugget, 0)
  Y_no_trend = simu$Y
  
  ## trend calculation
  
  trend = matrix(NA,K,N)
  for ( i in 1:N){
  trend[,i] = X[,1]*B[i]
  }
  
  
  ## scale trend in order to have random and fixed part in the same space
  trend = t(sapply(1:K,function(i){
    trend[i,] * diff(range(Y_no_trend[i,])) + mean(Y_no_trend[i,])
  }))
  
  ## sum the two parts
  Y =   trend + Y_no_trend
  ## used if the trend effect is by pieces 
  Y[,sX==0] = Y_no_trend[,sX==0]
  sX_T = sX == 1
  
  # scale the GP between 0 and 1
  Y = t(apply(Y,1,scale_01))
  
  
  return(Y)
}