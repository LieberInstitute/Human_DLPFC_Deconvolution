
library("spatialLIBD")
library("here")
library("sessioninfo")

## load snRNA-seq DLPFC data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## drop Ambiguous nuc
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

table(sce$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     10894      1940     24809     11067

broad_modeling_results <- registration_wrapper(
  sce,
  var_registration = "cellType_broad_hc",
  var_sample_id = "Sample",
  gene_ensembl = "gene_id",
  gene_name = "gene_name",
  pseudobulk_rds_file = here("processed-data", "sce","sce_broad_pseudobulk.rds")
)

save(broad_modeling_results, here("processed-data","06_marker_genes", "05_broad_registration", "broad_modeling_results.Rdata"))

# slurmjobs::job_single(name = "05_broad_registration", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 05_broad_registration.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

