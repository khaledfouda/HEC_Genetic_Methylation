simu_GaSP = function(D,
                     sites,
                     beta,
                     nugget,
                     tau_0,
                     seed = as.numeric(Sys.time()),
                     discrete = FALSE,
                     mispe = F,
                     sinuzo = F,
                     V_save = NULL,
                     trend = F) {
  
  N = length(sites)
  K = length(beta)
  set.seed(seed)
  A = matrix(runif(K*(K-D),0,0.1),nrow = K, ncol = K-D )
  #A = rustiefel(K,K-D)
  
  if(D > 0){
    if (discrete){
      if(D!=2){
        stop("in discrete case only D=2 allowed")
      }
      x = rep(c(1,0),ceiling(K/2))[1:K]
      X = matrix(0,ncol = 2, nrow= K)
      X[x==1,1] = 1
      X[x==0,2] = 1
    }else {
      X = matrix(runif(D*K),ncol = D, nrow= K)
      #X = rustiefel(K,D)
    }
    
    Hx = matrix(solve(t(X)%*%X,t(X)),nrow = length(X)/K)
    A = A-X%*%Hx%*%A
    if(mispe){
      #cat("ok")
      set.seed(seed)
      A = matrix(runif(K*(K-D),min(A),max(A)),nrow = K, ncol = K-D )
    }
    
    A = cbind(X,A)
    
  }
  
  
  

  

  #sapply(1:ncol(A), function(i) sum(X[,2]*A[,i]))
  
  
  if (length(V_save) < (K * N)) {
    V = matrix(NA, nrow = K, ncol = N)
    R00 = abs(outer(sites, sites, '-'))
    #for ( d in (D+1):K){
    for (d in 1:K) {
      R = matern_5_2_kernel(R00, beta = beta[d])
      R_tilde = R + nugget[d] * diag(N)
      set.seed(seed)
      #V[d, ] = rcpp_rmvnorm_stable(1, R, rep(0, N))
      V[d, ] = rcpp_rmvnorm_stable(1, R_tilde, rep(0, N))
    }
  } else {
    V = matrix(NA, nrow = K, ncol = N)
    R00 = abs(outer(sites, sites, '-'))
    if(D>0){
      for (d in 1:D) {
        R = exp_kernel(R00, beta = beta[d])
        R_tilde = R + nugget[d] * diag(N)
        set.seed(seed)
        #V[d, ] = rcpp_rmvnorm_stable(1, R, rep(0, N))
        V[d, ] = rcpp_rmvnorm_stable(1, R_tilde, rep(0, N))
      }
      
      
      V[-c(1:D),]= V_save[(D+1):K,]
    }else {
      V = V_save[1:K,]
    }

  }

  
  if(D>0){
    if(discrete){
      if(D!=2){
        stop("in discrete case only D=2 allowed")
      }
      
      if(sinuzo){
        V[1,] =  sin(2*pi*sites)
        V[2,] = -sin(2*pi*sites)
      }
    }else{
      if( trend ){
        print("trend ok")
        V[1,] = 3.5*sin(2*pi*sites)
      }
    }
    
  }
  
  

  
  # 
   Y_hat = A %*% V
  # 
  X = A[,-c((D+1):K)]
  A = A[,c((D+1):K)]
  B = V[-c((D+1):K),]
  V = V[c((D+1):K),]
  
  
  Y = matrix(rnorm(length(Y_hat),Y_hat,tau_0),nrow = nrow(Y_hat))
  
  mod_simu = list(
    Y = Y,
    sites = sites,
    A = A,
    V = V,
    tau_0 = tau_0,
    nugget = nugget, 
    beta = beta,
    B = B,
    X = X
  )
  return(mod_simu)
}