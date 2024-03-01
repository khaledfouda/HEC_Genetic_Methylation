get_dmr_regions_case1 <-
 function(p_values,
          sites,
          alpha = 1e-4,
          floor_by = 1e3,
          return_seq = TRUE) {
  #' if middle_point = FALSE then the sites are divided into equal regions of size floor_by and a region is considered
  #' dmr if it has at least min_freq of points below alpha/N.
  #' set return_seq = TRUE to return a list of all possible sites in the dmr regions
  #' if middle_point = TRUE  then it considers each point of p_value smaller than alpha/N and
  #' takes floor_by/2 of points before and after it
  #' if return_seq = FALSE: then a list of sites where the p_value is smaller than alpha/N is returned for middle_point=TRUE
  #' and  a list of starting positions of chosen regions for middle_point = FALSE
  
  N = length(sites)
  
  loci = sort(sites[which(p_values < (alpha / N))])
  if (!return_seq)
   return(loci)
  
  sequen = c()
  step_size = round(floor_by / 2)
  end_site = 1
  for (i in 1:length(loci)) {
   start_site = max(end_site, loci[i] - step_size)
   end_site = loci[i] + step_size
   sequen = c(sequen, start_site:end_site)
  }
  return(sequen)
  
 }
