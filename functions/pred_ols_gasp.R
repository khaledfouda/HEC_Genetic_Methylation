
pred_ols_gasp = function(obj){
  
  
  
  for (ind in 1:length(obj))
  {
    assign(names(obj)[ind], obj[[ind]])
  }
  
  K = nrow(Y_obs)
  D = ncol(A) 
  N = ncol(Y_obs)
  n_index = apply(is.na(Y_obs),2,sum)==0
  predict_all_dlm = sd_all_dlm = matrix(0,D, length(scale_sites))
  
  for(d in 1:D){
    #print(d)
    fgasp.model=fgasp(input,A_output_t[d,],have_noise = T)
    
    pred_fast=predict(param=log(c(beta_record[d],eta_record[d])),object=fgasp.model,
                      testing_input=scale_sites)
    
    #pred_fast=Kalman_smoother(log(c(beta_record[i_dlm],eta_record[i_dlm])), index_obs, t(A_output_t[i_dlm,]),delta_x_all,var_est_record[i_dlm])
    predict_all_dlm[d,]=pred_fast@mean
    pred_var = pred_fast@var
    pred_var[pred_var<0] = 0
    sd_all_dlm[d,]=sqrt(pred_var)
  }
  
  pred_all_record=A%*%(predict_all_dlm)
  
  # pred_all_record[,1]
  # (A%*%(predict_all_dlm))[,1]
  # (output_mean-Y_k)[,1]
 
  
  pred_testing_record_cond = var_testing_record_cond = matrix(NA,K,length(scale_sites))
  
  
  for(i in 1:length(scale_sites)){
    if(i%%10000==0){
      print(i)
    }
    testing_ppl_i = is.na(Y_obs[,i])
    
    if(sum(testing_ppl_i) == nrow(Y_obs)){
      pred_testing_record_cond[testing_ppl_i,i]=pred_all_record[testing_ppl_i,i]
      var_predict=Sigma_hat[testing_ppl_i,testing_ppl_i]
      var_testing_record_cond[testing_ppl_i,i]=diag(var_predict)
      
    }else if(sum(testing_ppl_i)>0){
      
      Sigma_hat=A%*%diag(sd_all_dlm[,i]^2)%*%t(A)
      test = try({inv_sig = solve(Sigma_hat[!testing_ppl_i,!testing_ppl_i])})
      if (class(test)[1] == "try-error") {
        inv_sig = solve(Sigma_hat[!testing_ppl_i,!testing_ppl_i] +
                          diag(rep(1e-6, sum(!testing_ppl_i))))
      }
      useful_block=Sigma_hat[testing_ppl_i,!testing_ppl_i]%*%inv_sig
      
      pred_testing_record_cond[testing_ppl_i,i]=pred_all_record[testing_ppl_i,i]+
        useful_block%*%(Y_obs[!testing_ppl_i,i]-rowMeans_t_output[!testing_ppl_i] -
                          pred_all_record[!testing_ppl_i,i])+ rowMeans_t_output[testing_ppl_i]
      
      # pred_testing_record_cond[testing_ppl_i,i]=pred_all_record[testing_ppl_i,i]+
      #   useful_block%*%(Y_obs[!testing_ppl_i,i] - 
      #                     pred_all_record[!testing_ppl_i,i])
      
      
      var_predict=Sigma_hat[testing_ppl_i,testing_ppl_i]-useful_block%*%Sigma_hat[!testing_ppl_i,testing_ppl_i]
      var_testing_record_cond[testing_ppl_i,i]=diag(var_predict)
      
    }
    
    if(sum(!testing_ppl_i)>0){
      pred_testing_record_cond[!testing_ppl_i,i]=Y_obs[!testing_ppl_i,i]#pred_all_record[!testing_ppl_i,i]
      var_testing_record_cond[!testing_ppl_i,i] = 0
    }
    
    
  }
  
  #Y_means = Y_obs-obj$rowMeans_t_output
  
  #i = 1
  #Y_means[,i]-pred_testing_record_cond[,i]
  #Y_k[,i]
  #plot(X,(Y_means - pred_testing_record_cond)[,i])
  
  #pred_testing_record_cond[,i] + X*B[i]
  #Y_means[,i]
  
  #Hx%*%Y_k[,i]
  #(Hx%*%(Y_means[,i]-pred_testing_record_cond[,i]))
  
  #B_est = apply(Y_means-pred_testing_record_cond,2,function(y) Hx[,!is.na(y)]%*%na.omit(y))
  
  
  #Y_pred = (X%*%t(B_est) + pred_testing_record_cond + rowMeans_t_output)
  #Y_pred = (X%*%t(B_est) + pred_testing_record_cond + rowMeans_t_output)
  Y_pred = pred_testing_record_cond
  
  return(Y_pred)
}