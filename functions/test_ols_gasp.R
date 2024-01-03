test_ols_gasp = function(Y_obs,sites,X,tol_eig = 0 ){
  
  if(sum(sites != sort(sites))>0){
    stop("sites must be sorted")
  }
  
  if(length(ncol(X))<1){
    nrow_X = length(X)
  }else{
    nrow_X = nrow(X)
  }
  
  
  K = nrow(Y_obs)
  D = K - (length(X)/K)
  n_index = apply(is.na(Y_obs),2,sum)==0
  scale_sites = scale_01(sites)
  input = scale_sites[n_index]
  
  if (length(input) < K) {
    warning(
      "ncol(Y_obs) without NA is lower than K, the SVD of I-Px(Y) will not be of rank K"
    )
  }
  
  if (K != nrow_X) {
    stop(
      "nrow(X) must be equal to nrow(Y_obs)"
    )
  }
  
  output = Y_obs[,n_index]
  C=(max(input)-min(input))/ncol(output)
  #rowMeans_t_output=rowMeans(output)
  rowMeans_t_output=rep(0,K)#rowMeans(Y_obs,na.rm = T)#
  
  output_mean = output - rowMeans_t_output
  
  
  
  # res_test = rowMeans(apply(output_mean,2,function(y){
  #   mod = lm(y~X-1)
  #   return(summary(mod)$coefficients[,4]<=0.001)
  # }))

  # res_test = rep(0.9,length(X)/K)
  # 
  # if(sum(res_test>0.5)<1){
  #   if(length(A)<1){
  #     print("ok")
  #     svd_output=svd(output-rowMeans_t_output)
  #     A=(svd_output$u%*%diag(svd_output$d)/sqrt(ncol(output)))
  #     A_output_t = t(svd_output$v)*sqrt(ncol(output))
  #     X = NULL
  #   }
  # }else {
    
    #X = X[,res_test>0.5]
    svd_X = svd(X)
    eig_X = svd_X$d
    select_eig = eig_X >= tol_eig
    eig_X[eig_X<=0] = 1e-30
    D_svd = diag(eig_X[select_eig],sum(select_eig),sum(select_eig))
    X_svd = svd_X$u[,1:sum(select_eig)]%*%D_svd
    
    inv_D_svd_2 = diag(1/(eig_X[select_eig]^2),sum(select_eig),sum(select_eig))
    B_svd = inv_D_svd_2%*%t(X_svd)%*%output_mean
    
    
    
    # Hx = matrix(solve(t(X)%*%X,t(X)),nrow = length(X)/K)
    # B = Hx%*%output_mean
    # 
    # plot(input,B[1,])
    
    Y_k = X_svd%*%B_svd#X%*%B#
    X = X_svd
    B = B_svd
    
    

    #### attention on n'a pas pris en compte le centrage des data !!!!!!!!!!
    
    
    #(Hx%*%Y_k)[1] == B[1]
    
    svd_A = svd(output_mean-Y_k)
    eig_A = svd_A$d
    select_eig = eig_A > tol_eig
    eig_A[eig_A<=0] = 1e-30
    D_svd = diag(eig_A[select_eig],sum(select_eig),sum(select_eig))
    # A = -(svd_A$u[,1:sum(select_eig)]%*%D_svd/sqrt(ncol(output)))[,1:min(D,sum(select_eig))]
    # A_output_t = t(-svd_A$v[,1:min(D,sum(select_eig))]*sqrt(ncol(output)))
    A = -(svd_A$u[,1:sum(select_eig)]%*%D_svd/sqrt(ncol(output)))[,1:sum(select_eig)]
    A_output_t = t(-svd_A$v[,1:sum(select_eig)]*sqrt(ncol(output)))
    
    
    A = cbind(X,A)
    A_output_t = rbind(B,A_output_t)
  #}
  
  obj = list(
    Y_obs = Y_obs,
    input = input,
    output = output,
    C = C,
    rowMeans_t_output = rowMeans_t_output,
    scale_sites = scale_sites,
    A = A,
    A_output_t = A_output_t,
    X = X,
    #Hx = Hx, 
    res_test= 0#res_test
      )
  
  return(obj)
}
