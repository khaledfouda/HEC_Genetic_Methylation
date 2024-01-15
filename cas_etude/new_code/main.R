
#----------------
# The following code reproduce our results
setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude")
source("new_code/load_files.R")
#------------------------------------------
chromosomes = paste0("chr", c(7,8,11,12,17))
alphas = c(.05, 0.1, 0.2)
#---------------------------------------------
# clean
transform_raw_to_feather(chromosome_to_retain = chromosomes)
#--------------------
condition = "DONOR_SEX == 'Male' & TISSUE_TYPE != 'Venous blood'"
condition = "is.na(DONOR_SEX) & TISSUE_TYPE != 'Venous blood'"
condition_post = "AGE != 2.5 & BONE_MARROW == 0"
note = "_MALE_BLOOD_"
chromosomes = c("chr17")
for(chr in chromosomes)
   combine_feathers_to_rds(chromosome = chr,condition_post = condition_post,note = note)
# OR for dmrseq
# transform_raw_to_feather2(chromosome_to_retain = chromosomes)
# for(chr in chromosomes)
#    combine_feathers_to_rds2(chromosome = chr)
#---------------------
# analyze
alphas = c(0.2)


# chromosome = chr; load_p_values = FALSE; alpha=alpha;
# min_freq=1; middle_point=TRUE; floor_by=1e3; plot=TRUE; note = note

for(chr in chromosomes){
   print(chr)
   for(alpha in alphas)
      compute_p_values_and_plot(chromosome = chr, load_p_values = T, alpha=alpha,
                                min_freq=1, middle_point=TRUE, floor_by=1e3, plot=TRUE, note = note)
}
# fit model: - warning: slow!
for(chr in chromosomes){
   for(alpha in alphas)
   res = run_model_on_chromosome(chr, subset=NA, min_freq = 1, no_cores = 5, middle_point=TRUE, floor_by=1e6,
                          alpha=alpha, Male.Only = TRUE, Age.Only = TRUE)
}
#------------------
# visualize results
results = get_results_data()
for(alpha in alphas)
   print(get_graph(results, alpha_val =  alpha))
#-----------
# END
