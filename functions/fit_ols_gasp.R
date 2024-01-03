fit_ols_gasp = function(obj){

  for (ind in 1:length(obj))
  {
    assign(names(obj)[ind], obj[[ind]])
  }

  D = ncol(A)

  beta_record=rep(0,D)
  eta_record=rep(0,D)
  val_record=rep(0,D)
  sigma2_record=rep(0,D)



  for(d in 1:D){
    fgasp.model=fgasp(input,A_output_t[d,], have_noise = T)

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


  obj = list(
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
    X = X,
   # Hx = Hx,
    res_test = res_test)

  return(obj)
}
