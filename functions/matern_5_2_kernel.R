matern_5_2_kernel<-function(d,beta){
  res=(1+sqrt(5)*beta*d + 5*beta^2*d^2/3 )*exp(-sqrt(5)*beta*d)
  res
}


