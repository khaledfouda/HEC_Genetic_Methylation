log_post<-function(param,object,C){
  log_lik_val=log_lik(param,object)
  beta=exp(param[1])
  eta=exp(param[2])  ####constraint
  a=1/2
  b=1
  logprior=approx_ref_matern_5_2(a,b,C,beta,  eta)
  log_lik_val+ logprior
}