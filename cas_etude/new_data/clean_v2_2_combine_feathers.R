
#---------------------------------------------------------------------------
combine_feathers_to_rds_2 <- function( chromosome = NULL,
                                    dat.info= readRDS("new_data/sample_info.rds"),
                                    cut_off = 1e6,
                                    save_all_chromosomes = FALSE
                                    ){
   # similar to combine_feathers_to_rds() but works for dmrseq.
   # saves X covariates, common positions, methylation matrix as a number and coverage matrix. 
   # the methylation value is  methylation/coverage
   
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
# example:
combine_feathers_to_rds_2(chromosome = "chr12") 
