approx_ref_matern_5_2<-function(a,b,C,beta_i,eta_i ){
  t=C*beta_i+eta_i  ###JR prior
  a*log(t)-b*t
}