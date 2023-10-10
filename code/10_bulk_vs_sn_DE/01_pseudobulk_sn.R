
library("SummarizedExperiment")
library("spatialLIBD")
library("here")
library("sessioninfo")

## load our DLPFC data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## drop Ambiguous nuc
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

table(sce$Sample)
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant 
# 3034        1889        2438        2708        5043        3200        2637        2701        2928        3182 
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid 
# 4492         944         921        4656        4009        3363        1927        4255        2120 

## revent issues with duplicated rownames
rownames(sce) <- rowData(sce)$gene_id


## Pseudobulk
sce_pb_sample <- registration_pseudobulk(
  sce,
  var_registration = "BrNum",
  var_sample_id = "Sample"
)

colData(sce_pb_sample)

colnames(sce_pb_sample) <- sce_pb_sample$Sample

## Save
save(sce_pb_sample, file = here("processed-data", "10_bulk_vs_sn_DE", "sce_pb_sample.Rdata"))

# slurmjobs::job_single(name = "01_pseudobulk_sn", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 01_pseudobulk_sn.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

