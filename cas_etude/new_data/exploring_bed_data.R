
#setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina")
setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude/")

# install rtracklayer to read bed files
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("rtracklayer")
#-----------

library(rtracklayer)
library(tidyverse)
library(readr)
library(magrittr)
library(skimr)
library(feather)
select <- dplyr::select

# test file
#wgbs_file1 <- './new_data/Gene_data/AML_BM_bp27_meth_methCpG.bed'
#df.test <- import.bed(wgbs_file1) %>% as.data.frame()

#table(df.test$seqnames) %>% as.data.frame() %>% arrange(desc(Var1)) 
#wgbs_file2 <- './new_data/Gene_data/TCD8_term_bp174_meth_methCpG.bed'

# take chr1 as a subset
#table(df$seqnames)



# Extracting methylation values
#df$site <- paste(df$seqnames, df$start, df$end, sep = "_")
#df$methylation <- sapply(strsplit(gsub("'", "", as.character(df$name)), "/"), 
#                         function(x) as.numeric(x[1]) / as.numeric(x[2]))
#head(df)

#view(head(df,10000))
#-------------------------------------------
transform_age <- function(age) {
   if (is.na(age) || age == "" || !grepl("^\\s*\\d+\\s*-\\s*\\d+\\s*$", age)) {
      return(NA)
   }
   
   age_range <- as.numeric(unlist(strsplit(age, "\\s*-\\s*")))
   print(age_range)
   return(mean(age_range))
}

chromosome_to_retain = paste0("chr", c(8:12,17))
gene_data_folder <- 'new_data/Gene_data/'
skip_check = TRUE
skip_sex = FALSE
save_all_chromosomes = TRUE
#---------------------------------------------------------------
run.the.slow.code = TRUE
#if(run.the.slow.code == TRUE){

transform_raw_to_feather <- function(chromosome_to_retain = NULL, 
                                     skip_sex = TRUE,
                                     skip_check = TRUE,
                                     gene_data_folder='new_data/Gene_data/',
                                     save_all_chromosomes = FALSE
                                     ){
         
   dat.info.full <- read_tsv("./new_data/samples.tsv", show_col_types = FALSE)
   
   dat.info.full %>%  
      select(sampleGroup, cellTypeShort, DONOR_ID, DONOR_AGE, DONOR_HEALTH_STATUS, DONOR_SEX,
             DISEASE, CELL_TYPE, TISSUE_TYPE,bedFile) %>%
      mutate(bedFile = paste0(gene_data_folder, bedFile)) %>%
      mutate(featherFile = sapply(bedFile, function(x) strsplit(x, "\\.")[[1]][1])) %>%
      mutate(across(-c(DONOR_AGE,bedFile, featherFile), as.factor)) %>%
      mutate(DONOR_AGE = sapply(DONOR_AGE, transform_age)) %>%
      filter(! is.na(DONOR_AGE)) ->
      dat.info  
   
   #skimr::skim(dat.info)
   #--------------------------------------------------------------
   XY.ratio <- function(seqnames){
      XY.count = as.character(seqnames[seqnames %in%  c("chrX", "chrY")]) %>% table()
      return(round(XY.count[1]/XY.count[2]))
   }
   #---------------------------------------------------------------------------------------
   dat.info$XY.ratio = NA
   for(i in 1:nrow(dat.info)){
      print(i)   
      gene.dat <- import.bed(dat.info$bedFile[i]) %>% as.data.frame()
      if(skip_sex == FALSE)
         dat.info$XY.ratio[i] = XY.ratio(gene.dat$seqnames)
      if(skip_check == FALSE)
         print(paste0( "end = start - 1? ", all(gene.dat %>% mutate(diff = (start - end)) %>% select(diff) == -1),
             ", # width != 2 is ", sum(gene.dat$width-2),
             ", X/Y = ", dat.info$XY.ratio[i]))
      
      if(save_all_chromosomes == TRUE){
         gene.dat %>% select(seqnames, start, score) %>% 
            write_feather(path=paste0(dat.info$featherFile[i],"_ALL_.feather"))
      }else{
         for(chromosome in chromosome_to_retain){
            gene.dat %>% filter(seqnames == chromosome) %>% select(start, score) %>% 
               write_feather(path=paste0(dat.info$featherFile[i],"_",chromosome,"_.feather"))
         }
      }
      print("...")
   }
   if(skip_sex == FALSE) 
      saveRDS(dat.info, "new_data/sample_info.rds")
   print("DONE.")
}
#-----------------------------------------------------------------------------
transform_raw_to_feather(save_all_chromosomes = TRUE)
#---------------------------------------------------------------------------
dat.info = readRDS("new_data/sample_info.rds")
run.the.second.slow.code = TRUE
skip_other_parts = TRUE
chromosome = "chr7"
cut_off = 1e6
for(chromosome in chromosome_to_retain){
   
#------------------------------------------------------------------------------
#if(run.the.second.slow.code == TRUE){

combine_feathers_to_rds <- function( chromosome = NULL,
                                    dat.info= readRDS("new_data/sample_info.rds"),
                                    skip_other_parts=TRUE,
                                    cut_off = 1e6,
                                    save_all_chromosomes = FALSE
                                    ){
   if(save_all_chromosomes == TRUE) chromosome = "ALL"
   
   X.dat <- 
      dat.info %>% 
      transmute(SAMPLE_ID = 1:nrow(dat.info),
                MALE = ifelse(is.na(DONOR_SEX), XY.ratio < 20, DONOR_SEX == "Male" ) + 0,
                AGE = DONOR_AGE,
                AML = (DISEASE %in% c("Acute myeloid leukemia", "Acute Myeloid Leukemia")) + 0,
                APL = (DISEASE == "Acute promyelocytic leukemia" ) + 0,
                BONE_MARROW = (TISSUE_TYPE == "Bone marrow") + 0,
                CH1_FILE = paste0(featherFile,"_",chromosome,"_.feather") ) 
   #-------------------------------------------------------------------------------
   # I need a list of common positions in Chromosome 1 for all patients.
   position_list = list()
   chromosome_list = list()
   
   X.dat$length = NA
   for(i in 1:nrow(X.dat)){
      feath.file = read_feather(X.dat$CH1_FILE[i])
      position_list[[i]] <- feath.file$start
      if(chromosome == "ALL")
         chromosome_list[[i]] <- position_list[[i]]
      
      position_list[[i]] <- read_feather(X.dat$CH1_FILE[i])$start
      X.dat$length[i] <- length(position_list[[i]])
   }
   
   position_list_2M = list()
   idx = 1
   for(i in 1:nrow(X.dat)){
      if(X.dat$length[i] > cut_off){
         position_list_2M[[idx]] = position_list[[i]]
         idx = idx + 1
      }
   }
   
   common_positions = Reduce(intersect, position_list_2M)
   
   lengths = unlist(lapply(position_list, length)) #%>% sort
   round(lengths / 1e3)
   round((lengths - mean(lengths))/ 1000  )
   #--------------------------------------------------------------------------
   # Y1 : Common positions only!
   N = length(common_positions)
   m = length(position_list_2M)
   print(N)
   X.dat_common = X.dat %>% mutate(length = (length - mean(length))/sd(length) )
   # My Y data will be a matrix 
   Y.dat_common <- matrix(NA, nrow = m, ncol = N)
   for(i in 1:nrow(X.dat_common)){
      (read_feather(X.dat_common$CH1_FILE[i]) %>%
         filter(start %in% common_positions) %>% 
         arrange(start))$score ->
         Y.dat_common[i,]
   }
   Y.dat_common = Y.dat_common / 1000
   saveRDS(sort(common_positions), paste0("new_data/sites_common_",chromosome,".rds"))
   saveRDS(as.matrix(select(X.dat_common,-CH1_FILE)), paste0("new_data/Xdat_common_",chromosome,".rds"))
   saveRDS(Y.dat_common, paste0("new_data/Ydat_common_",chromosome,".rds"))
   #---------------------------------------------------------------------------------
   if(skip_other_parts == FALSE){
      
      # Y2.1 :  Merge chromosomes together with different positions but ordered! Take the first non-NA sub-matrix
      # Include the two patients with short chromosome 1 (ie, length is less than 2M)
      lengths = unlist(lapply(position_list, length)) %>% sort
      N = min(lengths)
      m = length(position_list)
      X.dat_ordered1 = X.dat %>% mutate(length = (length - mean(length))/sd(length) )
      Y.dat_ordered1 = matrix(NA, nrow = m, ncol = N)
      for(i in 1:nrow(X.dat_ordered1)){
         (read_feather(X.dat_ordered1$CH1_FILE[i]) %>%
             arrange(start))$score[1:N] ->
            Y.dat_ordered1[i,]
      }
      Y.dat_ordered1 = Y.dat_ordered1 / 1000
      saveRDS(as.matrix(select(X.dat_ordered1,-CH1_FILE)), "new_data/Xdat_ordered1.rds")
      saveRDS(Y.dat_ordered1, "new_data/Ydat_ordered1.rds")
      #------------------------------------------------------------------------------------------------------
      # Y2.2 :  Merge chromosomes together with different positions but ordered! Take the first non-NA sub-matrix
      # consider chromosomes with more than 2M genes (ie, exclude two subjects)
      lengths = unlist(lapply(position_list_2M, length)) %>% sort
      N = min(lengths)
      m = length(position_list_2M)
      X.dat_ordered2 = X.dat_common %>% mutate(length = (length - mean(length))/sd(length) )
      Y.dat_ordered2 = matrix(NA, nrow = m, ncol = N)
      for(i in 1:nrow(X.dat_ordered2)){
         (read_feather(X.dat_ordered2$CH1_FILE[i]) %>%
             arrange(start))$score[1:N] ->
            Y.dat_ordered2[i,]
      }
      Y.dat_ordered2 = Y.dat_ordered2 / 1000
      saveRDS(as.matrix(select(X.dat_ordered2,-CH1_FILE)), "new_data/Xdat_ordered2.rds")
      saveRDS(Y.dat_ordered2, "new_data/Ydat_ordered2.rds")
   }
}
}
#------------------------------------------------------------------------------------------------------
# END
#######################
# gene.chr1$start %>% range()
# 
# #---------------------------------------------------------------------------------
# loaded.dat <- read_feather(dat.info$featherFile[i])
# 
# dat.info <- readRDS("new_data/sample_info.rds")
# 
# print(dat.info$DONOR_SEX[i])
# assign_sex(gen.dat.1$seqnames)
# 
# 
# 
# table((gen.dat.1 %>%
#    select(seqnames) %>%
#    filter(seqnames %in% c("chrX", "chrY")))$seqnames)
# 
