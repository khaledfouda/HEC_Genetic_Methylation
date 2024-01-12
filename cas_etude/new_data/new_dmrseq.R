setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")

#BiocManager::install("dmrseq")

library(BiocParallel)
library(tidyverse)
library(dmrseq)
chromosome = "chr7"
print(chromosome)
#for(chromosome in paste0("chr",c(8:12,17))){

chromosome = "chr7"
X = readRDS(paste0("new_data/Xdat_chr7_v2.rds")) %>% as.data.frame() #%>% 
Methylation = readRDS(paste0("new_data/Methylation_",chromosome,"_v2.rds")) %>% t()
Coverage = readRDS(paste0("new_data/Coverage_",chromosome,"_v2.rds")) %>% t()
sites = readRDS(paste0("new_data/sites_",chromosome,"_v2.rds"))
chr = rep(chromosome, nrow(Methylation))

for(chromosome in paste0("chr",c(8,11,12,17))){
   print(chromosome)
   Methylation = rbind(Methylation,t(readRDS(paste0("new_data/Methylation_",chromosome,"_v2.rds"))))
   Coverage = rbind(Coverage,t(readRDS(paste0("new_data/Coverage_",chromosome,"_v2.rds"))))
   chr = c(chr,rep(chromosome, nrow(Methylation)-length(sites) ))
   sites = c(sites,readRDS(paste0("new_data/sites_",chromosome,"_v2.rds")))
   
}

site_order = order(sites) 
Methylation = Methylation[site_order,]
Coverage = Coverage[site_order,]
sites = sites[site_order]
chr = chr[site_order]
table(chr)

female_indices = X$MALE==1
Methylation = Methylation[,female_indices]
Coverage = Coverage[, female_indices]
X = X[female_indices,]
# Methylation = readRDS(paste0("new_data/Methylation_",chromosome,"_v2.rds")) %>% t()
# Coverage = readRDS(paste0("new_data/Coverage_",chromosome,"_v2.rds")) %>% t()
# sites = readRDS(paste0("new_data/sites_",chromosome,"_v2.rds"))
# X = readRDS(paste0("new_data/Xdat_",chromosome,"_v2.rds")) %>% as.data.frame() #%>% 

dim(Methylation)
N = nrow(Methylation)
K = ncol(Methylation)
# chr = rep(chromosome, N)

if(length(chr) != nrow(Methylation) || length(sites) != nrow(Coverage)) {
   stop("Number of rows in Methylation and Coverage must match the length of chr and sites.")
}   
if(length(X$SAMPLE_ID) != ncol(Methylation) || length(X$SAMPLE_ID) != ncol(Coverage)) {
   stop("Number of columns in Methylation and Coverage must match the length of X$SAMPLE_ID.")
}
 

bs <- BSseq(chr = chr, pos = sites,
           M = Methylation, Cov = Coverage,
           sampleNames = X$SAMPLE_ID)    

pData(bs)$AGE <- X$AGE 
bs = sort(bs)
param <- SnowParam(workers = 1, type = "SOCK")
regions <- dmrseq(bs=bs, testCovariate="AGE", cutoff = .5,chrsPerChunk =1, BPPARAM = param)     
 
 # load example data
# data(BS.chr21)
# # the covariate of interest is the 'CellType' column of pData(BS.chr21)
# testCovariate <- 'CellType' 
# # run dmrseq on a subset of the chromosome (10K CpGs)
# regions <- dmrseq(bs=BS.chr21[240001:250000,],
#                   cutoff = 0.05,
#                   testCovariate=testCovariate)
# 
# 
# BiocManager::install(version = "3.18")
# BiocManager::install("bsseq")
# BiocManager::valid()
# BiocManager::version()
# 
# BiocManager::install("SummarizedExperiment")
# library(SummarizedExperiment)
# 
# BiocManager::install(c(
#    "bbmle", "bigalgebra", "brew", "brio", "cli", "cowplot", "curl", "dagitty", "data.table",
#    "datawizard", "DBI", "desc", "DT", "e1071", "emmeans", "fansi", "future", "future.apply",
#    "gdtools", "ggridges", "htmlwidgets", "httpuv", "igraph", "later", "lavaan", "maps", "markdown",
#    "matrixStats", "mixAK", "patchwork", "pkgbuild", "processx", "progress", "psych", "QuickJSR",
#    "ragg", "randomForestSRC", "Rcpp", "RCurl", "recipes", "s2", "sandwich", "sass", "segmented",
#    "stringi", "svglite", "tensorA", "tidygraph", "timeDate", "tseries", "vroom", "yaml"
# ), update = TRUE, ask = FALSE, force = TRUE)
# 
# infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
#                       package = 'bsseq')
# bismarkBSseq <- read.bismark(files = infile,
#                              rmZeroCov = TRUE,
#                              strandCollapse = FALSE,
#                              verbose = TRUE)
# bismarkBSseq@assays@data@listData[[3]][1:4,]
# 
# data(BS.chr21)
# 
# # reorder samples to create a null comparison 
# BS.null <- BS.chr21[1:20000,c(1,3,2,4)]
# 
# # add 100 DMRs
# BS.chr21.sim <- simDMRs(bs=BS.null, num.dmrs=100)
# 
# # bsseq object with original null + simulated DMRs
# show(BS.chr21.sim$bs)
# 
# 
# BS.chr21.sim$dmr.L[1:5]
