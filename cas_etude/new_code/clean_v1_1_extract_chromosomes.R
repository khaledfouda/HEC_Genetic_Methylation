
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

#-------------------------------------------
transform_age <- function(age) {
   if (is.na(age) || age == "" || !grepl("^\\s*\\d+\\s*-\\s*\\d+\\s*$", age)) {
      return(NA)
   }
   
   age_range <- as.numeric(unlist(strsplit(age, "\\s*-\\s*")))
   print(age_range)
   return(mean(age_range))
}
XY.ratio <- function(seqnames){
   XY.count = as.character(seqnames[seqnames %in%  c("chrX", "chrY")]) %>% table()
   return(round(XY.count[1]/XY.count[2]))
}
#---------------------------------------------------------------------------------------
transform_raw_to_feather <- function(chromosome_to_retain = NULL, 
                                     skip_sex = TRUE,
                                     skip_check = TRUE,
                                     gene_data_folder='new_data/Gene_data/',
                                     save_all_chromosomes = FALSE
                                     ){
   
   # Input :  a list of chromosome names to retain (each chromosome/subject pair is saved individually as a feather file).
   #        could also set save_all_chromosomes=TRUE to save all of the 26 chromosomes for each subject
   #        gene_data_folder is supposed to stay fixed
   #        skip_sex = FALSE: In order to predict the sex based on the numbers of X and Y chromosomes since some are missing.
   #        skip__check = FALSE: Verify that the width is always 2.
   
   
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
# using the function, or specify the chromosomes if needed
transform_raw_to_feather(save_all_chromosomes = TRUE)
#---------------------------------------------------------------------------
