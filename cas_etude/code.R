library(Jmisc)
library(DiceEval)
library(parallel)
library(plotly)
library(mgcv)
library(Rcpp)
library(RColorBrewer)
library(ggpubr)
library(ZIprop)
library(corrplot)


source("methyl_func.R")
sourceAll(path="../functions/")

## get data from fastgasp article
load("data/Methylation_level_dense_Chr1.Rda")
## get identified DMR area
load(file = "data/area_diff.Rdata")
load(file = "data/area_hypo.Rdata")



ind = 1:1e6
Y = t(seq.data.chr1[ind, ])
selec = c(1:2, 4:9, 12:14, 16:19)

Y = Y[selec, ]
colnames(Y) = rownames(Y) = NULL
X1 = c(rep("Gifford", 8), rep("roadmap1", 7))

sites = scale_01(Index_seq[ind])


N = ncol(Y)
K = nrow(Y)
methyl = Y
k_star = (1:K)[-c(1,9)]
Y_mean = apply(Y,2,mean)
X = fact2mat(X1)




#### set n_star ####

area_all = c(area_diff,area_hypo)

### 50 % of DMR area are taken
dmr_50_random = sort(unlist(lapply(area_all,function(x){
  sample(x,round(length(x)*0.5))
})))
dmr_50_random = which(ind %in% dmr_50_random )

### all DMR
ind_dmr = sort(unique(unlist(area_all)))
ind_dmr = which(ind %in% ind_dmr )



n_star_list = list(sort(unique(c(round(seq(1,N,length.out=round(N*0.9)))))), #scenario 1
                   sort(unique(c(round(seq(1,N,length.out=round(N*0.9))),dmr_50_random))), #scenario 2
                   sort(unique(c(round(seq(1,N,length.out=round(N*0.9))),ind_dmr))), #scenario 3
                   ind_dmr) #scenario 4

save(n_star, file = "data/n_star.Rdata")

sapply(n_star_list,length)/N


#### apply method ####
res = list()
length(res) = length(n_star_list)


for ( i in 1:4){
  n_star = n_star_list[[i]]
  if( i == 2){
    ind_na_sub = dmr_50_random
  }else if ( i %in% c(3,4)){
    ind_na_sub = ind_dmr
  }else if ( i == 1){
    ind_na_sub = NULL
  }
  res[[i]] = methyl_func(methyl, sites, k_star, n_star, X,ind_na_sub)
  save(res, file = "data/res_dmr_1e6_bis.Rdata")
}

res



