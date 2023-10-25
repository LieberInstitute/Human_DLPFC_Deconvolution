
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

# The variable to be tested must be a fixed effect
form <- ~library_type + (1|BrNum) + mitoRate + rRNA_rate + totalAssignedGene

#### RUN DE ####
source(here("code", "09_bulk_DE","run_DE.R"))

DREAM_library_type <- map2(rse_list, names(rse_list), function(rse, prep_name){
  
  message(Sys.time(), ' - Running DE ', prep_name)
  
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
  try(write.csv(outDE, file = here(data_dir, paste0("DREAM_library-type_gene_",prep_name,".csv")), row.names = TRUE))
  
  return(outDE)
})

save(DREAM_library_type, file = here(data_dir, "DREAM_library-type.Rdata"))

# slurmjobs::job_single(name = "08_DREAM_library-type", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 08_DREAM_library-type.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
