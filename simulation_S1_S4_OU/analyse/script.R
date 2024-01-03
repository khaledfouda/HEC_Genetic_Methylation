##----------------------------------------------------------------
##                          lien couleur                         -
##----------------------------------------------------------------
# https://www.color-hex.com/







##---------------------------------------------------------------
##                        Library loads                         -
##---------------------------------------------------------------

library(plotly)
library(DiceEval)
library(ggnewscale)
library(bannerCommenter)
library(data.table)
library(tidyr)
library(grDevices)
library(latex2exp)
library(Cairo)


op = par()


##---------------------------------------------------------------
##                          data load                           -
##---------------------------------------------------------------



source("boxplot_1.R")
source("scatter_1.R")
source("scatter_2.R")
source("scatter_3.R")
source("scatter_4.R")
source("barplot_1.R")


simu = c("discrete", "continuous" , "discrete_mispe", "discrete_mispe_block")


fig_time = T
save_pdf = T


for ( i in 1:4){
  
  
  filename = simu[i]
  load(file = paste(c("../data_res/res_",filename,"_OU.Rdata"),collapse = ""))
  
  
  ##---------------------------------------------------------------
  ##                            GRAPH 1                           -
  ##---------------------------------------------------------------
  name_fig = paste0("S",i)
  
  discrete = ifelse(i==2, F, T)
  
  
  
  l = 10
  
  x <- res[[l]]$X
  if(discrete){
    x <- res[[l]]$X[,1]
  }
  
  
  
  
  proc_a = which.min(x[11:length(x)])+10
  proc_b = which.max(x[11:length(x)])+10
  
  Y = res[[l]]$Y
  sites = res[[l]]$sites
  Y_obs = res[[l]]$Y_obs
  
  ind_na = is.na(Y_obs[50,])
  
  Y_GaSP = res[[l]]$Y_gasp
  Y_OLS_GaSP = res[[l]]$Y_ols_gasp
  
  
  
  n_rep = length(res)
  
  rmse = t(sapply(1:n_rep,function(i){
    ind_na = is.na(res[[i]]$Y_obs)
    Y = res[[i]]$Y
    sapply(res[[i]][10:11], function(al){
      if(length(al)==length(Y)){
        ind_pb = !is.na(al[ind_na])
        return(RMSE(Y[ind_na][ind_pb], al[ind_na][ind_pb]))
      }else{
        return(NA)
      }
    })
  }))
  
  
  colnames(rmse) = c("Y_lmcc","Y_gasp")
  
  
  
  boxplot_1(
    value = rmse,
    
    save_pdf = save_pdf,
    path = paste(c("figures/rmse_", name_fig, "_OU.eps"), collapse = ""),
    width = 10,
    height = 6,
    
    colorTitle = "black",
    sizeTitle = 25,
    formeTitle = "bold",
    colorAxe = "black",
    sizeAxe = 20,
    formeAxe = "bold.italic",
    textSize = 20,
    title = paste(c("Root Mean Square Error (S",i,")"),collapse = ""),
    
    col_boxplot = "black",
    fill_box = c("#E69F00","#45b9d2"),
    transparence_box = 0.9
  )
  
  
  
  
  
  
  
  ##---------------------------------------------------------------
  ##                            GRAPH 2                           -
  ##---------------------------------------------------------------
  
  
  if(discrete){
    scatter_4(
      sites = sites,
      Y = Y,
      Y_obs = Y_obs,#Y_OLS_GaSP,#,
      proc_a = proc_a,
      proc_b = proc_b,
      x = x,
      
      save_pdf = save_pdf,
      path = paste(c("figures/proc_obs_",name_fig,"_OU.eps"),collapse = ""),
      width = 10,
      height = 6,
      
      colorTitle = "black",
      sizeTitle = 25,
      formeTitle = "bold",
      
      colorAxe = "black",
      sizeAxe = 20,
      formeAxe = "bold",
      textSize = 20,
      Title = paste(c("Two samples (S",i,")"),collapse = ""),
      
      low = "#349be8",
      high = "#cc0000",
      point_size = 5
    )
  }else {
    scatter_3(
      sites = sites,
      Y = Y,
      Y_obs = Y_obs,
      proc_a = proc_a,
      proc_b = proc_b,
      x = x,
      
      save_pdf = save_pdf,
      path = paste(c("figures/proc_obs_",name_fig,"_OU.eps"),collapse = ""),
      width = 10,
      height = 6,
      
      colorTitle = "black",
      sizeTitle = 25,
      formeTitle = "bold",
      
      colorAxe = "black",
      sizeAxe = 20,
      formeAxe = "bold",
      
      textSize = 20,
      Title = paste(c("Two samples (S",i,")"),collapse = ""),
      
      low = "#349be8",
      high = "#cc0000",
      point_size = 5)
  }
}







##---------------------------------------------------------------
##                            GRAPH 3                           -
##---------------------------------------------------------------
#
# scatter_2(
#   sites = sites,
#   ind_na = ind_na,
#   Y = Y,
#   Y_GaSP = Y_GaSP,
#   Y_OLS_GaSP = Y_OLS_GaSP,
#   proc = proc_a,
#
#   save_pdf = F,
#   path = paste(c("figures/proc_res_2_",name_fig,".pdf"),collapse = ""),
#   width = 20,
#   height = 10,
#
#   colorTitle = "black",
#   sizeTitle = 20,
#   formeTitle = "bold",
#
#   colorAxe = "black",
#   sizeAxe = 15,
#   formeAxe = "bold",
#
#   textSize = 16,
#
#   col_gasp = "#f0d8c9",
#   col_ols_gasp = "#f46d61",
#   point_size = 2,
#   point_shape = 19
#
# )
#
#
#
# scatter_2(
#   sites = sites,
#   ind_na = ind_na,
#   Y = Y,
#   Y_GaSP = Y_GaSP,
#   Y_OLS_GaSP = Y_OLS_GaSP,
#   proc = proc_b,
#
#   save_pdf = F,
#   path = paste(c("figures/proc_res_2_",name_fig,".pdf"),collapse = ""),
#   width = 20,
#   height = 10,
#
#   colorTitle = "black",
#   sizeTitle = 20,
#   formeTitle = "bold",
#
#   colorAxe = "black",
#   sizeAxe = 15,
#   formeAxe = "bold",
#
#   textSize = 16,
#
#   col_gasp = "#f0d8c9",
#   col_ols_gasp = "#f46d61",
#   point_size = 2,
#   point_shape = 17
#
# )



##---------------------------------------------------------------
##                            TAbleau temps                     -
##---------------------------------------------------------------

if(fig_time){
  load(file = "../data_res/res_discrete_OU.Rdata")
  tab_time = as.vector(sapply(res, function(x) x$time[1:2]))
  load(file = "../data_res/res_continuous_OU.Rdata")
  tab_time = c(tab_time,as.vector(sapply(res, function(x) x$time[1:2])))
  load(file = "../data_res/res_discrete_mispe_OU.Rdata")
  tab_time = c(tab_time,as.vector(sapply(res, function(x) x$time[1:2])))
  load(file = "../data_res/res_discrete_mispe_block_OU.Rdata")
  tab_time = c(tab_time,as.vector(sapply(res, function(x) x$time[1:2])))
  tab_time = data.table(tab_time,
                        simulations = rep(paste0("S",1:4),each = 2*length(res)),
                        Method = rep(names(res[[1]])[c(10,11)],4*length(res)))
  names(tab_time)[1] = "Time"
  tab_time = tab_time[,mean(Time),by = list(simulations,Method)]
  names(tab_time)[3] = "Time"

  tab_time$Method = factor(tab_time$Method, levels = unique(tab_time$Method))
  levels(tab_time$Method) = c("Y_lmcc","Y_gasp")

  barplot_1(
    tab_time,
    save_pdf = T,
    path = paste(c("figures/time_OU.eps"), collapse = ""),
    width = 10,
    height = 4.5,

    colorTitle = "black",
    sizeTitle = 25,
    formeTitle = "bold",
    colorAxe = "black",
    sizeAxe = 20,
    formeAxe = "bold.italic",
    textSize = 20,
    title = "Average computation time (S1 to S4)",

    col_boxplot = "black",
    fill_box = c("#E69F00","#45b9d2"),
    transparence_box = 0.9
  )


}

toto = apply(rmse,2,mean)

cat("diff RMSE mean:",(toto[1]-toto[2])/toto[2])


