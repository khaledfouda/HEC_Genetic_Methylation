
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
BiocManager::install("ShortRead")
wgbs_file1 <- '../tmp/GSM1204463_BiSeq_cpgMethylation_BioSam_1500_HepG2_304072.BiSeq.bed'
df.test <- import.bed(wgbs_file1) %>% as.data.frame()

table(df.test$seqnames) %>% as.data.frame() %>% arrange(desc(Var1)) 
filter(df.test, seqnames == "chr7")$start -> tmp1
filter(df.test, seqnames == "chr8")$start -> tmp2
filter(df.test, seqnames == "chr1")$start -> tmp3
par(mfrow=c(3,1))
hist(tmp1)
hist(tmp2)
inters <- Reduce(intersect, list(tmp1, tmp2))
hist(inters)
#wgbs_file2 <- './new_data/Gene_data/TCD8_term_bp174_meth_methCpG.bed'
length(inters)
length(tmp1)
# take chr1 as a subset
#table(df$seqnames)
range(tmp1)
sort(tmp3)[1:20]
Reduce(intersect, list(tmp3, Index_seq)) -> tmp4
length(Index_seq) - length(tmp4)

filter(df.test, start == 96)
min(df.test$start)

df.test$width[5]
df.test[1:5,]
as.integer(strsplit(strsplit(df.test$name, "'")[[1]][2],"/")[[1]])

cov.meth = strsplit(trimws(df.test$name,whitespace = '\''),'/') 
Meth = sapply(cov.meth, function(x) as.integer(x[1])) 
Cov =  sapply(cov.meth, function(x) as.integer(x[2]))

Meth <- mutate(df.test, Meth = as.integer(strsplit(trimws(name,whitespace = '\''),'/')[[1]][1]))

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
XY.ratio <- function(seqnames){
   XY.count = as.character(seqnames[seqnames %in%  c("chrX", "chrY")]) %>% table()
   return(round(XY.count[1]/XY.count[2]))
}
#---------------------------------------------------------------------------------------
transform_raw_to_feather_2 <- function(chromosome_to_retain = NULL, 
                                     skip_sex = TRUE,
                                     skip_check = TRUE,
                                     gene_data_folder='new_data/Gene_data/',
                                     save_all_chromosomes = FALSE
                                     ){
         
   dat.info <-
      read_tsv("./new_data/samples.tsv", show_col_types = FALSE) %>% 
      select(sampleGroup, cellTypeShort, DONOR_ID, DONOR_AGE, DONOR_HEALTH_STATUS, DONOR_SEX,
             DISEASE, CELL_TYPE, TISSUE_TYPE,bedFile) %>%
      mutate(bedFile = paste0(gene_data_folder, bedFile)) %>%
      mutate(featherFile = sapply(bedFile, function(x) strsplit(x, "\\.")[[1]][1])) %>%
      mutate(across(-c(DONOR_AGE,bedFile, featherFile), as.factor)) %>%
      mutate(DONOR_AGE = sapply(DONOR_AGE, transform_age)) %>%
      mutate(XY.ratio = NA) %>% 
      filter(! is.na(DONOR_AGE))
   
   #skimr::skim(dat.info)
   #--------------------------------------------------------------
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
         print("isn't implemented yet")
         return()
         gene.dat %>% select(seqnames, start, score) %>% 
            write_feather(path=paste0(dat.info$featherFile[i],"_ALL_.feather"))
      }else{
         for(chromosome in chromosome_to_retain){
            gene.dat %>% filter(seqnames == chromosome) %>% select(start, name) %>% 
               write_feather(path=paste0(dat.info$featherFile[i],"_",chromosome,"_v2.feather"))
         }
      }
      print("...")
   }
   if(skip_sex == FALSE) 
      saveRDS(dat.info, "new_data/sample_info.rds")
   print("DONE.")
}
#-----------------------------------------------------------------------------
chromosome_to_retain = paste0("chr", c(7,8,11,12,17))
transform_raw_to_feather_2(chromosome_to_retain)
#---------------------------------------------------------------------------
combine_feathers_to_rds_2 <- function( chromosome = NULL,
                                    dat.info= readRDS("new_data/sample_info.rds"),
                                    skip_other_parts=TRUE,
                                    cut_off = 1e6,
                                    save_all_chromosomes = FALSE
                                    ){
   if(save_all_chromosomes == TRUE){
      print("isn't implemented yet")
      return()
      chromosome = "ALL"
   }   
   
   X.dat <- 
      dat.info %>% 
      transmute(SAMPLE_ID = 1:nrow(dat.info),
                MALE = ifelse(is.na(DONOR_SEX), XY.ratio < 20, DONOR_SEX == "Male" ) + 0,
                AGE = DONOR_AGE,
                AML = (DISEASE %in% c("Acute myeloid leukemia", "Acute Myeloid Leukemia")) + 0,
                APL = (DISEASE == "Acute promyelocytic leukemia" ) + 0,
                BONE_MARROW = (TISSUE_TYPE == "Bone marrow") + 0,
                CH1_FILE = paste0(featherFile,"_",chromosome,"_v2.feather") ) 
   #-------------------------------------------------------------------------------
   # I need a list of common positions in Chromosome 1 for all patients.
   position_list = list()
   chromosome_list = list()
   
   X.dat$length = NA
   for(i in 1:nrow(X.dat)){
      if(chromosome == "ALL")
      {
         feath.file = read_feather(X.dat$CH1_FILE[i])
         position_list[[i]] <- feath.file$start
         chromosome_list[[i]] <- feath.file$seqnames
      }else
         position_list[[i]] <- read_feather(X.dat$CH1_FILE[i])$start
      
      X.dat$length[i] <- length(position_list[[i]])
      
   }
   
   if(chromosome != "ALL"){
      
      position_list_2M = list()
      idx = 1
      
      for(i in 1:nrow(X.dat)){
         if(X.dat$length[i] > cut_off){
            position_list_2M[[idx]] = position_list[[i]]
            idx = idx + 1
         }
      }
   }else{
      position_list_2M = position_list
   }
   position_list = NULL   
   common_positions = Reduce(intersect, position_list_2M)
   
   
   #lengths = unlist(lapply(position_list, length)) #%>% sort
   #round(lengths / 1e3)
   #round((lengths - mean(lengths))/ 1000  )
   #--------------------------------------------------------------------------
   # Y1 : Common positions only!
   N = length(common_positions)
   m = length(position_list_2M)
   print(N)
   X.dat = X.dat %>% mutate(length = (length - mean(length))/sd(length) )
   # My Y data will be a matrix 
   
   Coverage <- matrix(NA, nrow=m, ncol = N)
   Methylation <- matrix(NA, nrow=m, ncol = N)
   for(i in 1:nrow(X.dat)){
      meth.file <-   
         read_feather(X.dat$CH1_FILE[i]) %>%
         filter(start %in% common_positions) %>%
         arrange(start)
      cov.meth = strsplit(trimws(meth.file$name,whitespace = '\''),'/') 
      Methylation[i,] = as.integer(sapply(cov.meth, function(x) x[1])) 
      Coverage[i,] =  as.integer(sapply(cov.meth, function(x) x[2]))
   }
   saveRDS(sort(common_positions), paste0("new_data/sites_",chromosome,"_v2.rds"))
   saveRDS(as.matrix(select(X.dat,-CH1_FILE)), paste0("new_data/Xdat_",chromosome,"_v2.rds"))
   saveRDS(Methylation, paste0("new_data/Methylation_",chromosome,"_v2.rds"))
   saveRDS(Coverage, paste0("new_data/Coverage_",chromosome,"_v2.rds"))
   if(chromosome == "ALL")
      saveRDS(order(common_positions), paste0("new_data/chromosomes_common_",chromosome,".rds"))
   #---------------------------------------------------------------------------------
}
#------------------------------------------------------------------------------------------------------
chromosome = "chr7"
chromosome_to_retain = paste0("chr", c(7,8,11,12,17))
for(chromosome in chromosome_to_retain)
   combine_feathers_to_rds_2(chromosome) 
#------------------------------------------------------------------------------

#-------------------------------------------------------------
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
