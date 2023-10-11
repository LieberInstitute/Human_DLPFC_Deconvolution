
library("SummarizedExperiment")
library("edgeR")
library("variancePartition")
library("purrr")
library("here")
library("jaffelab")
library("sessioninfo")

#### Set up ####
## dirs
# plot_dir <- here("plots", "09_bulk_DE", "08_DREAM_library-type")
# if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## dirs
data_dir <- here("processed-data", "09_bulk_DE", "09_DREAM_library-prep")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

rse_list <- map(splitit(rse_gene$library_type), ~rse_gene[,.x])
map_int(rse_list, ncol)
# polyA RiboZeroGold 
# 55           55 

map(rse_list, ~table(.x$library_prep))
# $polyA
# 
# Bulk Cyto  Nuc 
# 19   18   18 
# 
# $RiboZeroGold
# 
# Bulk Cyto  Nuc 
# 19   19   17

## pairwise comaprison of library_preps
prep <- c("Bulk", "Cyto", "Nuc")
prep_pairs <- map(prep, ~prep[prep != .x])
names(prep_pairs) <- map_chr(prep_pairs, ~paste0(.x, collapse = "_"))

# The variable to be tested must be a fixed effect
form <- ~library_prep + (1|BrNum) + mitoRate + rRNA_rate + totalAssignedGene

#### RUN DE ####
DREAM_library_prep <- map2(rse_list, names(rse_list), function(rse, type_name){
  
  map2(prep_pairs, names(prep_pairs), function(pair, pair_name){
    message(Sys.time(), ' - Running DE ', type_name, "+", pair_name)
    
    rse <- rse[,rse$library_prep %in% pair]
    message("Subsample to:", ncol(rse))
    
    dge = DGEList(assays(rse)$counts)
    dge = calcNormFactors(dge)
    
    metadata = as.data.frame(colData(rse))
    
    param = SnowParam(4, "SOCK", progressbar=TRUE)
    # estimate weights using linear mixed model of dream
    message(Sys.time(), ' - voomWithDreamWeights')
    vobjDream = voomWithDreamWeights(dge, form, metadata, BPPARAM=param)
    
    # Fit the dream model on each gene
    # By default, uses the Satterthwaite approximation for the hypothesis test
    message(Sys.time(), ' - Dream')
    fitmm = dream(vobjDream, form, metadata)
    message(Sys.time(), ' - eBayes')
    fitmm = eBayes(fitmm)
    
    # Examine design matrix
    head(fitmm$design, 3)
    
    outDE <- topTable( fitmm, coef="library_typeRiboZeroGold", number=Inf , sort.by = "none")
    
    message(Sys.time(), " - Saving")
    try(write.csv(outDE, file = here(data_dir, paste0("DREAM_library-prep_gene_",type_name,"_",pair_name,".csv")), row.names = FALSE))
    
    return(outDE)
  })
  
})

save(DREAM_library_type, file = here(data_dir, "DREAM_library-type.Rdata"))

# slurmjobs::job_single(name = "09_DREAM_library-prep", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 09_DREAM_library-prep.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
