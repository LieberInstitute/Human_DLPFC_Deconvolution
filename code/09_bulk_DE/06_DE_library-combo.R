
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("purrr")
library("here")
library("jaffelab")
library("sessioninfo")

#### Set up ####

## plot dir
plot_dir <- here("processed-data", "09_bulk_DE", "06_DE_library-combo")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## data dir
data_dir <- here("processed-data", "09_bulk_DE", "06_DE_library-combo")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#### load data ####
features <- c("gene", "exon", "jx", "tx")

rse_paths <- here("processed-data","rse", paste0("rse_", features, ".Rdata"))
names(rse_paths) <- features

## for testing
# rse_paths <- rse_paths["gene"]

rse_list <- lapply(rse_paths, function(x) get(load(x)))

#### build model ####
pd <- as.data.frame(colData(rse_list$gene))
pd$Sample <- factor(pd$Sample)
mod <- model.matrix(~Sample + library_type + library_prep + mitoRate + rRNA_rate + totalAssignedGene, data = pd)
dim(mod)
# [1] 110  22

colnames(mod)

#### RUN DE ####
source(here("code", "09_bulk_DE","run_DE.R"))

DE_library_combo <- map2(rse_list, names(rse_list), function(rse, feat_name){
    message(Sys.time(), ' - Running DE ', feat_name)
    
    outDE <- run_DE(rse = rse,
                    model = mod, 
                    coef = c("library_typeRiboZeroGold", "library_prepCyto", "library_prepNuc"), 
                    run_voom = feat_name != "tx",
                    save_eBayes = TRUE,
                    plot_name = here(plot_dir, paste0("DE_library-combo_", feat_name, "_",prep_name,".pdf"))
                      )
    
    message(Sys.time(), " - Saving")
    try(write.csv(outDE$topTable, file = here(data_dir, paste0("DE_library-combo_", feat_name, ".csv")), row.names = FALSE))
    
    return(outDE)
})


save(DE_library_combo, file = here(data_dir, "DE_library-combo.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
