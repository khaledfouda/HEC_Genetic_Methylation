
pred_fgasp = function(obj){



  for (ind in 1:length(obj))
  {
    assign(names(obj)[ind], obj[[ind]])
  }


  K = nrow(Y_obs)
  D = ncol(A)

  ncol_X = 0
  if(length(X)>0){
    ncol_X = 1:(length(X)/K)
  }

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
      inv_sig = Sigma_hat[!testing_ppl_i,!testing_ppl_i]

      test = try({inv_sig = solve(Sigma_hat[!testing_ppl_i,!testing_ppl_i])},silent = T)
      if (class(test)[1] == "try-error") {
        inv_sig = solve(Sigma_hat[!testing_ppl_i,!testing_ppl_i] +
                          diag(rep(1e-6, sum(!testing_ppl_i))))
      }
      useful_block=Sigma_hat[testing_ppl_i,!testing_ppl_i]%*%inv_sig

      pred_testing_record_cond[testing_ppl_i,i]=pred_all_record[testing_ppl_i,i]+
        useful_block%*%(Y_obs[!testing_ppl_i,i]-rowMeans_t_output[!testing_ppl_i] -
                          pred_all_record[!testing_ppl_i,i])+ rowMeans_t_output[testing_ppl_i]


      var_predict=Sigma_hat[testing_ppl_i,testing_ppl_i]-useful_block%*%Sigma_hat[!testing_ppl_i,testing_ppl_i]
      var_testing_record_cond[testing_ppl_i,i]=diag(var_predict)

    }

    if(sum(!testing_ppl_i)>0){
      pred_testing_record_cond[!testing_ppl_i,i]=Y_obs[!testing_ppl_i,i]
      var_testing_record_cond[!testing_ppl_i,i] = 0
    }


  }


  Y_pred = pred_testing_record_cond

  return(Y_pred)
}
