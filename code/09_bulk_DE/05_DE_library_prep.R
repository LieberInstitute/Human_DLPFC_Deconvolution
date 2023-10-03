
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("purrr")
library("here")
library("jaffelab")
library("sessioninfo")

#### Set up ####

## dirs
data_dir <- here("processed-data", "09_bulk_DE", "04_DE_library-type")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
features <- c("gene", "exon", "jx", "tx")

rse_paths <- here("processed-data","rse", paste0("rse_", features, ".Rdata"))
names(rse_paths) <- features

all(file.exists(rse_paths))

rse_list <- lapply(rse_paths, function(x) get(load(x)))

## edit rse list: factor sample, split library_prep

rse_list <- map(rse_list, function(rse){
  rse <- map(splitit(rse$library_prep), ~rse[,.x])
  return(rse)
})

#### build model ####
mod <- map(rse_list$gene, function(rse){
  pd <- as.data.frame(colData(rse))
  pd$Sample <- factor(pd$Sample)
  mod <- model.matrix(~Sample + library_type, data = pd)
  return(mod)
})

map(mod, dim)
# $Bulk
# [1] 38 20
# 
# $Cyto
# [1] 37 20
# 
# $Nuc
# [1] 35 20

colnames(mod[[1]])

#### RUN DE ####
source(here("code", "09_bulk_DE","run_DE.R"))

DE_out <- map2(rse_list, names(rse_list), function(rse, feat_name){
  pmap(list(rse_prep = rse, prep_name = names(rse), mod = mod), function(rse_prep, prep_name, mod){
    
    message(Sys.time(), ' - Running DE ', feat_name, " + ", prep_name)
    stopifnot(ncol(rse_prep) == nrow(mod))
    
    outDE <- run_DE(rse = rse_prep, model = mod, coef = "library_typeRiboZeroGold", run_voom = feat_name != "tx")
    
    message(Sys.time(), " - Saving")
    write.csv(outDE, file = here(data_dir, paste0("DE_library-type_", feat_name, "_",prep_name,".csv")), row.names = FALSE)
    
    return(outDE)
  })
})

head(DE_out$gene$Bulk)

save(DE_out, file = here(data_dir, "DE_library-type.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
