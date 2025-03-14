
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("purrr")
library("here")
library("jaffelab")
library("sessioninfo")

#### Set up ####

## plot dir
plot_dir <- here("plots", "09_bulk_DE", "05_DE_library-prep")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## data dirs
data_dir <- here("processed-data", "09_bulk_DE", "05_DE_library-prep")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
features <- c("gene", "exon", "jx", "tx")

rse_paths <- here("processed-data","rse", paste0("rse_", features, ".Rdata"))
names(rse_paths) <- features

all(file.exists(rse_paths))

## for testing
# rse_paths <- rse_paths["gene"]

rse_list <- lapply(rse_paths, function(x) get(load(x)))

## edit rse list: factor sample, split library_prep
rse_list <- map(rse_list, function(rse){
  rse <- map(splitit(rse$library_type), ~rse[,.x])
  return(rse)
})

#### build model ####
mod <- map(rse_list$gene, function(rse){
  pd <- as.data.frame(colData(rse))
  pd$Sample <- factor(pd$Sample)
  mod <- model.matrix(~Sample + library_prep + mitoRate + rRNA_rate + totalAssignedGene, data = pd)
  return(mod)
})

map(mod, dim)
# $polyA
# [1] 55 21
# 
# $RiboZeroGold
# [1] 55 21

colnames(mod[[1]])

#### RUN DE ####
source(here("code", "09_bulk_DE","run_DE.R"))

DE_library_prep <- map2(rse_list, names(rse_list), function(rse, feat_name){
  pmap(list(rse_prep = rse, prep_name = names(rse), mod = mod), function(rse_prep, prep_name, mod){
    
    message(Sys.time(), ' - Running DE ', feat_name, " + ", prep_name)
    stopifnot(ncol(rse_prep) == nrow(mod))
    
    if(feat_name == "gene"){
      plot_name <- here(plot_dir, paste0("DE_library-prep_", feat_name, "_",prep_name,".pdf"))
    } else {
      plot_name <- NULL
    }
    
    outDE <- run_DE(rse = rse_prep, 
                    model = mod, 
                    coef = c("library_prepCyto", "library_prepNuc"), 
                    run_voom = feat_name != "tx", 
                    save_eBayes = TRUE,
                    plot_name = plot_name
    )
    
    message(Sys.time(), " - Saving")
    try(write.csv(outDE$topTable, file = here(data_dir, paste0("DE_library-type_", feat_name, "_",prep_name,".csv")), row.names = FALSE))
    
    return(outDE)
  })
})


save(DE_library_prep, file = here(data_dir, "DE_library-type.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
