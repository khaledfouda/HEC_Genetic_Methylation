
#----------------
# The following code reproduce our results
setwd("/mnt/campus/math/research/kfouda/main/HEC/Melina/latest/cas_etude")
source("new_code/load_files.R")
#------------------------------------------
chromosomes = paste0("chr", c(1:12,17))
note = "_subset_Blood"
alphas = c(1e-4)
alpha =alphas[1]
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
# P -values: 2 options
for(chr in chromosomes){
   print(chr)
   for(alpha in alphas) 
      compute_p_values_and_plot(chromosome = chr, load_p_values = F, alpha=alpha, Male.Only = FALSE,
                                min_freq=1, middle_point=TRUE, floor_by=1e3, plot=TRUE, note = note,
                                correction = correction) 
}
#----------
# Multivariate p-values
for(chr in chromosomes){
   print(chr)
   for(alpha in alphas) 
      compute_p_values_and_plot_multivar(chromosome = chr, load_p_values = T,
                                         alpha=alpha, Male.Only = FALSE,
                                min_freq=1, middle_point=TRUE, floor_by=1e3, plot=T, note = note,
                                correction = correction) 
}

#-------------
# combine generated images to a single pdf

image_folder <- "./graphs/"
image_files1 <- list.files(image_folder, pattern = "manhattan_plot_Multivar_1_.*\\.png$", full.names = TRUE)
image_files2 <- list.files(image_folder, pattern = "manhattan_plot_Multivar_2_.*\\.png$", full.names = TRUE)


pdf("./graphs/combined_p_values_Multivar.pdf", width = 8.5, height = 11)

for(i in 1:length(image_files1)) {
   plot.new()
   par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
   plot(0:1, 0:1, type = "n", xlab = "", ylab = "", axes = FALSE, ann = FALSE)
   
   x_left = (8.5 - 8) / 2 / 8.5
   x_right = 1 - x_left
   y_bottom = 0.5 + (0.5 - 8 / 22) / 2
   y_top = 1 - (0.5 - 8 / 22) / 2
   img <- png::readPNG(image_files1[i])
   rasterImage(img, x_left, y_bottom, x_right, y_top)
   
   y_bottom = (0.5 - 8 / 22) / 2
   y_top = 0.5 - (0.5 - 8 / 22) / 2
   img <- png::readPNG(image_files2[i])
   rasterImage(img, x_left, y_bottom, x_right, y_top)
}
dev.off()



#---------------------------------------------------
# fit model: - warning: slow!
for(chr in chromosomes){
   print(chr)
   for(alpha in alphas)
   res = run_model_on_chromosome(chr, subset=NA, min_freq = 1, no_cores = 2, 
                                 middle_point=TRUE, floor_by=1e3,
                          alpha=alpha, Male.Only = FALSE, Age.Only = FALSE,
                          note = note, correction=correction) 
}    
#-----------------------------------------------------------------------
# fit model: - warning: slow!
for(chr in paste0("chr",c(1:4,6:12,17))){
#   chr = "chr5"
   print(chr)
      res = run_model_on_chromosome_Multivar(chr, subset=NA, min_freq = 1, no_cores = 2, 
                                    middle_point=TRUE, floor_by=1e3, 
                                    alpha=alpha, Male.Only = FALSE, Age.Only = FALSE, 
                                    note = note, correction=correction) 
}    



 
#----------------------------------------------------------------------------
# fit model but only on scenario 4
for(chr in chromosomes){
   print(chr)
   for(alpha in alphas)
      
   res = run_model_on_chromosome_dmr(chr, Age.Only=TRUE, Male.Only=FALSE,alpha=alpha,
                                        min_freq=1, min_k=2,subset=NA, no_cores=3,
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
  