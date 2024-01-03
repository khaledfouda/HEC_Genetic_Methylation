fit_fgasp = function(Y_obs,sites,A = NULL, D = nrow(Y_obs)){
  
  if(sum(sites != sort(sites))>0){
    stop("sites must be sorted")
  }
  
  K = nrow(Y_obs)
  D = K
  n_index = apply(is.na(Y_obs),2,sum)==0
  scale_sites = scale_01(sites)
  input = scale_sites[n_index]
  
  output = Y_obs[,n_index]
  C=(max(input)-min(input))/ncol(output)
  rowMeans_t_output=rowMeans(Y_obs,na.rm = T)#rep(0,K)#
  
  if(D>ncol(output)){
    stop("D must be less than or equal to N - length(n_star) ")
  }
  
  if(length(A)<1){
    svd_output=svd(output-rowMeans_t_output)
    A=(svd_output$u%*%diag(svd_output$d)/sqrt(ncol(output)))
    A = A[1:K,1:D]
    A_output_t = t(svd_output$v)*sqrt(ncol(output))
  }
  
  # test = try({inv_A = solve(t(A)%*%A)})
  # if(class(test)[1] == "try-error"){
  #   inv_A = solve(t(A)%*%A+diag(1e-6,ncol(A),ncol(A)))
  # }
  # A_output_t=inv_A%*%t(A)%*%(output-rowMeans_t_output)
  # 
  
  D = ncol(A)
  
  
  beta_record=rep(0,D)
  eta_record=rep(0,D)
  val_record=rep(0,D)
  sigma2_record=rep(0,D)
  
  
  
  
  for(d in 1:D){
    #print(d)
    fgasp.model=fgasp(input,A_output_t[d,],have_noise = T)
    
    test = try({
      tt_all <-
        optim(c(log(1 / C), 1), function(par, object)
          return(log_post(par, object, C)), object = fgasp.model, lower = c(-100,-15),
          method = "L-BFGS-B",
          control = list(fnscale = -1, maxit = 30))
    })
    
    count = 0
    while (class(test)[1] == "try-error" & count<=100) {
      count = count + 5
      test = try({
        tt_all <-
          optim(c(log(1 / C), 1), function(par, object)
            return(log_post(par, object, C)), object = fgasp.model, lower = c(-100/count,-15),
            method = "L-BFGS-B",
            control = list(fnscale = -1, maxit = 30))
      })
    }

    #tt_all
    val_record[d]=tt_all$value
    
    beta_record[d]=exp(tt_all$par)[1]
    eta_record[d]=exp(tt_all$par)[2]
    
    sigma2_record[d]=Get_log_det_S2(param=log(c(beta_record[d],eta_record[d])),fgasp.model@have_noise,fgasp.model@delta_x,
                                    fgasp.model@output,fgasp.model@kernel_type)[[2]]/fgasp.model@num_obs
    
  }
    
  
obj_fgasp = list(
  Y_obs = Y_obs,
  input = input,
  output = output,
  rowMeans_t_output = rowMeans_t_output,
  scale_sites = scale_sites,
  A = A,
  A_output_t = A_output_t,
  beta_record = beta_record,
  eta_record = eta_record,
  val_record = val_record,
  sigma2_record = sigma2_record,
  X = NULL)

  return(obj_fgasp)
}