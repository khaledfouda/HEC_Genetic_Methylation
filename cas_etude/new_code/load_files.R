setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")

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
library(tidyverse)
library(magrittr)

library(doParallel)
library(foreach)

library(rtracklayer)
library(tidyverse)
library(readr)
library(magrittr)
library(gridExtra)
library(skimr)
library(feather)
select <- dplyr::select




sourceAll(path="../functions/")
select <- dplyr::select
source("new_code/clean_v1_1_extract_chromosomes.R")
source("new_code/clean_v1_2_combine_feathers.R")
source("new_code/clean_v2_1_extract_chromosomes.R")
source("new_code/clean_v2_2_combine_feathers.R")

source("new_code/analyze_1_compute_p_values.R")
source("new_code/analyze_2_compute_dmr.R")
source("new_code/analyze_3_fit_models.R")
source("new_code/analyze_4_fit_models_scenario4.R")
#source("new_code/results_1_analyze.R")