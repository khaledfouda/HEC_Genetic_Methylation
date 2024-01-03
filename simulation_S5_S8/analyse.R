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
library(tidyr)
library(Cairo)

op = par()

##---------------------------------------------------------------
##                          data load                           -
##---------------------------------------------------------------



source("../simulation_S1_S4/analyse/boxplot_1.R")
source("../simulation_S1_S4/analyse/barplot_1.R")

fig_time = T
save_pdf = T
S = c(c(5:8),"5_OU")
for ( i in 1:5){
  
  file_name = paste(c("data/res_mcci",i,".Rdata"), collapse = "")
  load(file = file_name)
  
  
  
  
  
  n_rep = max(which(sapply(res,length)>0))
  
  rmse = t(sapply(1:n_rep,function(i){
    ind_na = is.na(res[[i]]$Y_obs)
    Y = res[[i]]$Y
    sapply(res[[i]][c("Y_ols_gasp","Y_mcci")], function(al){
      if(length(al)==length(Y)){
        ind_pb = !is.na(al[ind_na])
        return(RMSE(Y[ind_na][ind_pb], al[ind_na][ind_pb]))
      }else{
        return(NA)
      }
    })
  }))
  
  colnames(rmse) = c("Y_lmcc","Y_mcci")
  
  boxplot_1(
    value = rmse,
    
    save_pdf = save_pdf,
    path = paste0(c("figures/rmse_mcci_",i,".eps"),collapse =""),
    width = 10,
    height = 6,
    
    colorTitle = "black",
    sizeTitle = 25,
    formeTitle = "bold",
    colorAxe = "black",
    sizeAxe = 20,
    formeAxe = "bold.italic",
    textSize = 20,
    title = paste(c("Root Mean Square Error (S",S[i],")"),collapse = ""),
    
    col_boxplot = "black",
    fill_box = c("#E69F00","#c9a0dc"),
    transparence_box = 0.9
  )
}



# tab_time = as.vector(sapply(res, function(x) x$time[1:2]))
# tab_time = data.table(tab_time,
#                       Method = rep(names(res[[1]][c("Y_ols_gasp","Y_mcci")]),length(res)),
#                       simulations = "")
# names(tab_time)[1] = "Time"
# tab_time = tab_time[,sum(Time),by = list(simulations,Method)]
# names(tab_time)[3] = "Time"
#
# tab_time$Method = factor(tab_time$Method, levels = unique(tab_time$Method))
#
# barplot_1(
#   tab_time,
#   save_pdf = save_pdf,
#   path = paste0(c("figures/time_mcci_",i,".pdf"),collapse =""),
#   width = 10,
#   height = 6,
#
#   colorTitle = "black",
#   sizeTitle = 25,
#   formeTitle = "bold",
#   colorAxe = "black",
#   sizeAxe = 20,
#   formeAxe = "bold.italic",
#   textSize = 20,
#
#   col_boxplot = "black",
#   fill_box = c("#E69F00","#c9a0dc"),
#   transparence_box = 0.9
# )




if(fig_time){
  tab_time_vec = NULL
  vec_n_rep = rep(NA,4)
  for(i in 1:4){
    print(i)
    file_name = paste(c("data/res_mcci",i,".Rdata"), collapse = "")
    load(file = file_name)
    vec_n_rep[i] =n_rep = max(which(sapply(res,length)>0))
    tab_time_vec = c(tab_time_vec,as.vector(sapply(res[1:n_rep], function(x) x$time[1:2])))
  }


  tab_time = data.table(tab_time_vec/60,
                        simulations = rep(paste0("S",1:4),2*vec_n_rep),
                        Method = as.factor(rep(c("Y_ols_gasp","Y_mcci"),sum(vec_n_rep))))
  names(tab_time)[1] = "Time"
  tab_time = tab_time[,mean(Time),by = list(simulations,Method)]
  names(tab_time)[3] = "Time"

  tab_time$Method = factor(tab_time$Method, levels = unique(tab_time$Method))
  levels(tab_time$Method) = c("Y_lmcc","Y_mcci")


  barplot_1(
    tab_time,
    save_pdf = save_pdf,
    path = paste(c("figures/time_mcci.eps"), collapse = ""),
    width = 10,
    height = 4.5,

    colorTitle = "black",
    sizeTitle = 25,
    formeTitle = "bold",
    colorAxe = "black",
    sizeAxe = 20,
    formeAxe = "bold.italic",
    textSize = 20,
    ylab = "Time in minutes",
    title = "Average computation time (S5 to S8)",

    col_boxplot = "black",
    fill_box = c("#E69F00","#45b9d2"),
    transparence_box = 0.9
  )


}
