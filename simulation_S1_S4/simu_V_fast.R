library(mvmeta)
library(Jmisc)
library(Rcpp)

sourceAll("functions/")

sourceAll("../functions/")

N = 1000
K = 100
n_rep = 100

t1 = tic()
V_list = list()
length(V_list) = n_rep
for ( seed in 1:n_rep){
  print(seed)
  set.seed(seed)
  sites = sort(runif(N))
  set.seed(seed)
  beta = runif(K,10,1000)
  set.seed(seed)
  nugget = runif(K,0.001,0.05)
  
  V = matrix(NA, nrow = K, ncol = N)
  R00=abs(outer(sites,sites,'-'))
  #for ( d in (D+1):K){
  for ( d in 1:K){
    R=matern_5_2_kernel(R00,beta=beta[d])
    R_tilde=R+nugget[d]*diag(N)
    set.seed(seed)
    #V[d,] = rcpp_rmvnorm_stable(1,R,rep(0,N))
    V[d, ] = rcpp_rmvnorm_stable(1, R_tilde, rep(0, N))
  }
  V_list[[seed]] = V
}

save(V_list, file = "data_res/V_list_100.Rdata")
