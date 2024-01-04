setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")


library(DiceEval)
library(knitr)
library(kableExtra)
library(data.table)




load(file = "data/res_dmr_1e6.Rdata")
tabl = rbindlist(lapply(res,function(toto){
  data.table(toto)
}))
tabl$time = round(tabl$time/60,2)

tabl = data.frame(tabl)
tabl = round(tabl,3)
tabl

#rownames(tabl) = paste0(c("OLS_GaSP","GaSP","Null"),rep(1:3,each = 3))
rownames(tabl) = paste0(c("OLS_GaSP","GaSP","Null"),rep(1:4,each = 3))
#colnames(tabl) = c("RMSE", "R2", "Computation time in sec")
colnames(tabl) = c("RMSE", "R2", "RMSE_DMR", "R2_DMR","Computation time in sec")

kbl(tabl, booktabs = T, format = "latex",row.names = T) %>%
  pack_rows("NA on the whole domain (900000)", 1, 3)%>%
  pack_rows("NA on a part of the specific areas and on the rest of the domain (902000)", 4, 6)%>%
  pack_rows("NA on specific areas and on the rest of the domain (904000)", 7, 9)%>%
  pack_rows("NA only on specific areas (4000)", 10, 12)










tabl
# the following file does not exist.
data = fread("data/descr_sample_methyl.csv")
tabl = data[,c(2,3,5,6)]
kbl(tabl, booktabs = T, format = "latex",row.names = T)



