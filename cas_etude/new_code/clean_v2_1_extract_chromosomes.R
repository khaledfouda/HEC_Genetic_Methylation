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
   #' The difference between this and transform_raw_to_feather_2 is that this extracts the mythylation is "4/10" instead of a value.
   #' This is needed to use the dmrseq package.
   #' Input: chromosome_to_retain: A list of chromosomes to retain from each participant, each chromosome is saved
   #'          to its own file as "sample_name_chromosome_v2.feather     
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
# Example use:
# chromosome_to_retain = paste0("chr", c(7,8,11,12,17))
# transform_raw_to_feather_2(chromosome_to_retain = chromosome_to_retain)
#------------------------------------------------------------------------------------