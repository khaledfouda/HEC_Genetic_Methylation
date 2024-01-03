#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/MacDoc/PostDocHecUqam/code_article_lm/simulations")

library(mvmeta)
library(Jmisc)
library(Rcpp)

## functions to create GP
sourceAll("functions/")

## functions for gls 
#sourceCpp(file="../src/beta_gls.cpp")

## functions to estimate GP
sourceAll("../functions/")


## data to accelerate simu 
load(file="data_res/V_list_100_OU.Rdata")

n_rep = 100



N = 1000
K = 100
k_star = c(11:K)

method = 1:2





#### S1 ####

sX = rep(1,N)

#
res = list()
length(res) = n_rep
for ( seed in 1:n_rep){
  print(seed)
  set.seed(seed)
  n_star = sort(sample(1:N,round(N*0.9)))
  #n_star = c(26:475,526:975)
  res[[seed]]  = run_simu(N,K,2,n_star,k_star,method ,sX,seed,discrete = T,sinuzo = T,
                          V = V_list[[seed]])
}

save(res, file = "data_res/res_discrete_OU.Rdata")
# # 

# # 

#### S2 ####
res = list()
length(res) = n_rep
for ( seed in 1:n_rep){
  print(seed)
  set.seed(seed)
  n_star = sort(sample(1:N,round(N*0.9)))
  #n_star = c(26:475,526:975)
  res[[seed]]  = run_simu(N,K,1, n_star, k_star, method ,sX,seed,discrete=F,
                          V = V_list[[seed]], trend = T)
}

save(res, file = "data_res/res_continuous_OU.Rdata")
# #
# # #
# #

#### S3 ####

res = list()
length(res) = n_rep
for ( seed in 1:n_rep){
  print(seed)
  set.seed(seed)
  n_star = sort(sample(1:N,round(N*0.9)))
  #n_star = c(26:475,526:975)
  res[[seed]]  = run_simu(N,K,2, n_star, k_star, method ,sX,seed,discrete=T, mispe=T, sinuzo = T,
                          V = V_list[[seed]])
}

save(res, file = "data_res/res_discrete_mispe_OU.Rdata")


#### S4 ####

res = list()
length(res) = n_rep
for ( seed in 1:n_rep){
  print(seed)
  set.seed(seed)
  #n_star = sort(sample(1:N,round(N*0.5)))
  n_star = c(26:475,526:975)
  #n_star = c(3:48,53:98)
  res[[seed]]  = run_simu(N,K,2, n_star, k_star, method ,sX,seed,discrete=T, sinuzo = T,mispe = T)
}

save(res, file = "data_res/res_discrete_mispe_block_OU.Rdata")


#### other test not in the article ####

# #
# #
# res = list()
# length(res) = n_rep
# for ( seed in 1:n_rep){
#   print(seed)
#   set.seed(seed)
#   n_star = sort(sample(1:N,round(N*0.9)))
#   #n_star = c(26:475,526:975)
#   res[[seed]]  = run_simu(N,K,1, n_star, k_star, method ,sX,seed,discrete=F, mispe=T, sinuzo = F)#,
#                           #V = V_list[[seed]])
# }
# 
# save(res, file = "data_res/res_continuous_mispe.Rdata")
# #
# #
# res = list()
# length(res) = n_rep
# for ( seed in 1:n_rep){
#   print(seed)
#   set.seed(seed)
#   n_star = sort(sample(1:N,round(N*0.9)))
#   #n_star = c(26:475,526:975)
#   res[[seed]]  = run_simu(N,K,2, n_star, k_star, method ,sX,seed,discrete=T,wrong_X = T, sinuzo = T,
#                           V = V_list[[seed]])
# }
# 
# save(res, file = "data_res/res_discrete_wrong_X.Rdata")
# 
# 
# res = list()
# length(res) = n_rep
# for ( seed in 1:n_rep){
#   print(seed)
#   set.seed(seed)
#   n_star = sort(sample(1:N,round(N*0.9)))
#   #n_star = c(26:475,526:975)
#   res[[seed]]  = run_simu(N,K,1, n_star, k_star, method ,sX,seed,discrete=F,wrong_X = T, sinuzo = F,
#                           V = V_list[[seed]])
# }
# 
# save(res, file = "data_res/res_continuous_wrong_X.Rdata")
# 
# 
# # 
