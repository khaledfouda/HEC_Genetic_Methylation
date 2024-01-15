combine_feathers_to_rds <- function( chromosome = NULL,
                                    dat.info= readRDS("new_data/sample_info.rds"),
                                    cut_off = 1e6,
                                    condition_pre = NA,
                                    condition_post = NA,
                                    save_all_chromosomes = FALSE,
                                    note = ""
                                    ){
   # Warning: If save_all_chromosomes is chosen then all chromosomes will be saved to the same matrix!!
   # run the function on each chromosome individually.
   # Saves a Y matrix of mythylation values, X matrix of covariates, and common sites vector for the positions
   # it only considers common sites and individuals with number of genes per chromosome larger than the cut_off
   
   if(save_all_chromosomes == TRUE) chromosome = "ALL"
   
   if(!is.na(condition_pre))
      dat.info <- dat.info %>% filter(eval(parse(text=condition_pre)))
   
   X.dat <- 
      dat.info %>% 
      transmute(SAMPLE_ID = 1:nrow(dat.info),
                MALE = ifelse(is.na(DONOR_SEX), XY.ratio < 20, DONOR_SEX == "Male" ) + 0,
                AGE = DONOR_AGE,
                AML = (DISEASE %in% c("Acute myeloid leukemia", "Acute Myeloid Leukemia")) + 0,
                APL = (DISEASE == "Acute promyelocytic leukemia" ) + 0,
                BONE_MARROW = (TISSUE_TYPE == "Bone marrow") + 0,
                CH1_FILE = paste0(featherFile,"_",chromosome,"_.feather") )
   
   if(!is.na(condition_post))
      X.dat <- X.dat %>% filter(eval(parse(text=condition_post)))
   #-------------------------------------------------------------------------------
   # I need a list of common positions in Chromosome 1 for all patients.
   position_list = list()
   chromosome_list = list()
   
   X.dat$length = NA
   for(i in 1:nrow(X.dat)){
      feath.file = read_feather(X.dat$CH1_FILE[i])
      position_list[[i]] <- feath.file$start
      if(chromosome == "ALL")
         chromosome_list[[i]] <- feath.file$seqnames
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
   saveRDS(sort(common_positions), paste0("new_data/sites_common_",chromosome,note,".rds"))
   saveRDS(as.matrix(select(X.dat_common,-CH1_FILE)), paste0("new_data/Xdat_common_",chromosome,note,".rds"))
   saveRDS(Y.dat_common, paste0("new_data/Ydat_common_",chromosome,note,".rds"))
   if(chromosome == "ALL")
      saveRDS(order(common_positions), paste0("new_data/chromosomes_common_",chromosome,".rds"))
   #---------------------------------------------------------------------------------
}
#---------------------------------------------------------------
# Example:
#combine_feathers_to_rds(chromosome = "chr12")
