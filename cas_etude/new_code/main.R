
#----------------
# The following code reproduce our results
setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude")
source("new_code/load_files.R")
#------------------------------------------
chromosomes = paste0("chr", c(6:1))#c(7:12,17))
note = "_subset_Blood"
alphas = c(1e-4)
correction =  function(alpha,N) alpha
#chromosomes = c("chr1")
#alphas = c(.05, 0.1, 0.2)
#---------------------------------------------
# clean
#transform_raw_to_feather(chromosome_to_retain = chromosomes)
#--------------------
condition_pre = "DONOR_SEX == 'Male' & TISSUE_TYPE != 'Venous blood'"
condition_pre = "is.na(DONOR_SEX) & TISSUE_TYPE != 'Venous blood'"

#(BONE_MARROW == 1)
condition_post = "(AGE != 2.5)  & (BONE_MARROW == 0) & (! SAMPLE_ID %in% c(7,17,56,36))"

#note = "_subset_BONE_"



for(chr in chromosomes)
   combine_feathers_to_rds(chromosome = chr,condition_post = condition_post,note = note)
# OR for dmrseq
# transform_raw_to_feather2(chromosome_to_retain = chromosomes)
# for(chr in chromosomes)
#    combine_feathers_to_rds2(chromosome = chr)
#---------------------
# analyze


# chromosome = chr; load_p_values = T; alpha=alpha;
# min_freq=1; middle_point=TRUE; floor_by=1e3; plot=TRUE; note = note

for(chr in chromosomes){
   print(chr)
   for(alpha in alphas) 
      compute_p_values_and_plot(chromosome = chr, load_p_values = F, alpha=alpha, Male.Only = FALSE,
                                min_freq=1, middle_point=TRUE, floor_by=1e3, plot=TRUE, note = note,
                                correction = correction) 
}
# fit model: - warning: slow!
for(chr in chromosomes){
   print(chr)
   for(alpha in alphas)
   res = run_model_on_chromosome(chr, subset=NA, min_freq = 1, no_cores = 2, 
                                 middle_point=TRUE, floor_by=1e3,
                          alpha=alpha, Male.Only = FALSE, Age.Only = FALSE,
                          note = note, correction=correction) 
}    

# fit model but only on scenario 4
for(chr in chromosomes){
   print(chr)
   for(alpha in alphas)
      
   res = run_model_on_chromosome_dmr(chr, Age.Only=TRUE, Male.Only=FALSE,alpha=alpha,
                                        min_freq=1, min_k=2,subset=NA, no_cores=2,
                                        floor_by=1e3, middle_point=TRUE, note = note,
                                        correction = correction)
}

#------------------
# visualize results
results = get_results_data()  
for(alpha in alphas)    
   print(get_graph(results, alpha_val =  alpha))    
#----------- 
# END
# BiocManager::install("GenomicFeatures")
# library(GenomicFeatures) 
# BiocManager::install("RMariaDB")
# library(rtracklayer)
# # Example: Download human gene annotations 
# txdb <- makeTxDbFromEnsembl(organism = "Homo sapiens", release = 75) 
#  
# 
# genes <- genes(txdb)
# head(genes)
# export.bed(genes, "gene_annotations.bed")
# library(GenomicRanges) 
# wgbs_data <- import.bed("path_to_your_wgbs_file.bed") 
# # Find overlaps
# overlaps <- findOverlaps(wgbs_data, genes) 
# 
# # Extract overlapping WGBS sites
# overlapping_wgbs <- wgbs_data[queryHits(overlaps), ]
# 
# # Extract overlapping genes
# overlapping_genes <- genes[subjectHits(overlaps), ]
# library(Gviz)
# 
# # Example visualization
# plotTracks(list(GenomeAxisTrack(), AnnotationTrack(genes), DataTrack(wgbs_data)), 
#            from = <start_region>, to = <end_region>)
  