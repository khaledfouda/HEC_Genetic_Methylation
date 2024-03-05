library(tidyverse)

note = ""
correction <- function(alpha, N)
 alpha / N
alpha = .05
region_length <- 1e3
dmr.info <- methyl.info <- data.frame()

# load data
source("methyl_func.R")
sourceAll(path = "../functions/")

## get data from fastgasp article
load("data/Methylation_level_dense_Chr1.Rda")
## get identified DMR area
load(file = "data/area_diff.Rdata")
load(file = "data/area_hypo.Rdata")
Y = t(seq.data.chr1)
selec = c(1:2, 4:9, 12:14, 16:19)
Y = Y[selec,]
colnames(Y) = rownames(Y) = NULL
X1 = c(rep("Gifford", 8), rep("roadmap1", 7))
X = c(rep(1, 8), rep(0, 7))
sites = Index_seq
p_values = readRDS(paste0("new_data/case1_p_values",note, ".rds"))
#------------

N = length(sites)
alpha = correction(alpha, N)
loci = sort(sites[which(p_values < alpha)])
sequen = numeric(length(loci))
step_size = round(region_length / 2)
end_site = 1

for (i in 1:length(loci)) {
 start_site = max(end_site, loci[i] - step_size)
 end_site = loci[i] + step_size
 sequen[i] = length(which(sites %in%  (start_site:end_site)))
}

data.frame(site = loci,
           region_length = sequen) ->
 dmr.info

#-------------------------
# Get Methylation info
dmr_regions <- get_dmr_regions_case1(p_values,
                                     sites,
                                     alpha,
                                     floor_by = region_length,
                                     return_seq = TRUE)

ind_dmr = sort(which(sites %in% dmr_regions))

Y = colMeans(Y)
# data.frame(chromosome = chromosome, site = sites, Methylation = Y[ind_dmr]) %>%
data.frame(
 site = sites,
 Methylation = Y,
 dmr = 1:N
) %>%
 mutate(dmr = ifelse(dmr %in% ind_dmr, TRUE, FALSE)) ->
 methyl.info


saveRDS(dmr.info, paste0("new_data/case1_dmr_info_", note, ".rds"))
saveRDS(methyl.info, paste0("new_data/case1_methyl_info_", note, ".rds"))
