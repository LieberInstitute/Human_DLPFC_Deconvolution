
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
data_dir <- here("processed-data", "09_bulk_DE", "08_DREAM_library-type")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
load(here("processed-data","rse", "rse_gene.Rdata"), verbose = TRUE)

## edit rse list: factor sample, split library_prep

rse_list <- map(splitit(rse_gene$library_prep), ~rse_gene[,.x])

map(rse_list, dim)

#### build model ####
mod <- map(rse_list, function(rse){
  pd <- as.data.frame(colData(rse))
  pd$Sample <- factor(pd$Sample)
  mod <- model.matrix(~library_type + BrNum + mitoRate + rRNA_rate + totalAssignedGene, data = pd)
  # mod <- model.matrix(~library_type + Sample, data = pd) # Simple Model
  return(mod)
})

map(mod, dim)
# $Bulk
# [1] 38 14
# 
# $Cyto
# [1] 37 14
# 
# $Nuc
# [1] 35 14

colnames(mod[[1]])

# The variable to be tested must be a fixed effect
form <- ~library_type + (1|BrNum) + mitoRate + rRNA_rate + totalAssignedGene

#### RUN DE ####
source(here("code", "09_bulk_DE","run_DE.R"))

DREAM_library_type <- pmap(list(rse_prep = rse, prep_name = names(rse)), function(rse_prep, prep_name, mod){
  
  message(Sys.time(), ' - Running DE ', prep_name)
  
  dge = DGEList(assays(rse)$counts)
  dge = calcNormFactors( dge )
  
  metadata = as.data.frame(colData(rse))
  
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights(dge, form, metadata, BPPARAM=param )
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  fitmm = dream( vobjDream, form, metadata )
  fitmm = eBayes(fitmm)
  
  # Examine design matrix
  head(fitmm$design, 3)
  
  outDE <- topTable( fitmm, coef='Disease1', number=Inf , sort.by = "none")
  
  message(Sys.time(), " - Saving")
  try(write.csv(outDE, file = here(data_dir, paste0("DE_library-type_", feat_name, "_",prep_name,".csv")), row.names = FALSE))
  
  return(outDE)
})


save(DREAM_library_type, file = here(data_dir, "DREAM_library-type.Rdata"))

# slurmjobs::job_single(name = "01_pseudobulk_sn", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 01_pseudobulk_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
