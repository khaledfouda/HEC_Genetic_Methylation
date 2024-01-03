gls_gasp = function(formula , Y_obs, sites, X,D = nrow(Y_obs),sX_T = rep(TRUE,ncol(Y_obs))){
  
  K = nrow(Y_obs)
  N = ncol(Y_obs)
  
  rowMeans_obs=rowMeans(Y_obs,na.rm = T)
  
  Y_obs = Y_obs - rowMeans_obs
  
  X_vec = data.frame(sapply(1:ncol(X), function(i) rep(X[,i],sum(sX_T))))
  names(X_vec) = names(X)
  cols = names(which(sapply(X,class) == "factor"))
  X_vec[cols] <- lapply(X_vec[cols], factor)
  
  
  
  
  Y_gam = data.table(Y = as.vector(Y_obs[,sX_T]), 
                     sites = rep(sites[sX_T],each = K),
                     X_vec)
  
  mod_gam = gam(formula, data = Y_gam)
  Y_pred_gam = matrix(0,nrow=K,ncol=N)
  Y_pred_gam[,sX_T] = matrix(predict(mod_gam,Y_gam[,-1]),nrow=K)
  Y_no_trend = Y_obs -  Y_pred_gam
  obj_gls =  fit_fgasp(Y_no_trend ,sites, D = D)
  
  Y_no_trend = pred_fgasp(obj_gls)
  
  Xv = predict.gam(mod_gam,Y_gam[,-1],type="lpmatrix")

  A = obj_gls$A
  beta = obj_gls$beta_record
  eta = obj_gls$eta_record
  sigma2 = obj_gls$sigma2_record
  input = sites
  
  Y_hat = Y_gam$Y
  obs_index = apply(Y_obs[,sX_T],2,function(x) which(!is.na(x))-1)
  N_obs = sum(sapply(obs_index,length))
  ind_na = is.na(Y_hat)
  Xv_na = Xv[!ind_na,]
  

  
  beta_gls = coef_gls_cpp(Y_hat[!ind_na],Xv_na,obs_index, N_obs ,sites[sX_T],beta,eta,sigma2,A)
  

  Y_pred_gam_gls = matrix(0,nrow=K,ncol=N)
  Y_pred_gam_gls[,sX_T] = matrix(Xv%*%beta_gls,nrow=K)
  Y_GLS_GaSP = pred_fgasp(obj_gls) + Y_pred_gam_gls
  
  return(Y_GLS_GaSP + rowMeans_obs)
  
}




