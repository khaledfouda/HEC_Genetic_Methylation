library(tidyverse)

chromosome = "chr2"
note = "_subset_Blood"
correction <- function(alpha,N) alpha
alpha = 1e-4
region_length <- 1e3
dmr.info <- methyl.info <- data.frame()

run_code = FALSE
if(run_code){
 for(chromosome in paste0("chr",c(1:12,17))){
  print(chromosome)
  sites = readRDS(paste0("new_data/sites_common_",chromosome, note, ".rds"))
  p_values = readRDS(paste0("new_data/p_values_",chromosome, note, ".rds"))
  N = length(sites)
  alpha = correction(alpha, N)
  loci = sort(sites[which(p_values < alpha)])
  sequen = numeric(length(loci))
  step_size = round(region_length/2)
  end_site = 1
  
  for(i in 1:length(loci)){
   start_site = max( end_site, loci[i]-step_size)
   end_site = loci[i] + step_size
   sequen[i] = length(which(sites %in%  (start_site: end_site)))
  }
  
  data.frame(chromosome = chromosome, site=loci, region_length=sequen) %>%
   rbind(dmr.info) ->
   dmr.info
  
  #-------------------------
  # Get Methylation info
  dmr_regions = get_dmr_regions(p_values, sites, alpha, floor_by = region_length,
                                middle_point = T, return_seq = T)
  ind_dmr = sort(which(sites %in% dmr_regions))
  #sites = sites[ind_dmr]
  Y = readRDS(paste0("new_data/Ydat_common_",chromosome, note, ".rds"))
  Y = colMeans(Y)
  # data.frame(chromosome = chromosome, site = sites, Methylation = Y[ind_dmr]) %>%
  data.frame(chromosome = chromosome, site = sites, Methylation = Y, dmr = 1:N) %>%
   mutate(dmr = ifelse(dmr %in% ind_dmr, TRUE, FALSE)) %>% 
   rbind(methyl.info) ->
   methyl.info
  
 }
 saveRDS(dmr.info, paste0("new_data/dmr_info_",note, ".rds"))
 saveRDS(methyl.info, paste0("new_data/methyl_info_",note, ".rds"))