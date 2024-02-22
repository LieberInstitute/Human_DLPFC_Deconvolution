library("SingleCellExperiment")
library("here")
library("sessioninfo")

## Access Tran et al. 
# Download and save a local cache of the data available at:
# https://github.com/LieberInstitute/10xPilot_snRNAseq-human#processed-data
bfc <- BiocFileCache::BiocFileCache()
url <- paste0(
  "https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/",
  "SCE_DLPFC-n3_tran-etal.rda"
)
local_data <- BiocFileCache::bfcrpath(url, x = bfc)
# load sce.dlpfc.tran
load(local_data, verbose = TRUE)

sce.dlpfc.tran$cellType_broad <- gsub("_[A-Z]","",sce.dlpfc.tran$cellType)
sce.dlpfc.tran$cellType_broad <- gsub("Mural", "EndoMural", sce.dlpfc.tran$cellType_broad)

exp_cell_types <- c("Astro", "EndoMural", "Micro", "Oligo","OPC", "Excit", "Inhib")
sce.dlpfc.tran <- sce.dlpfc.tran[sce.dlpfc.tran$cellType_broad %in% exp_cell_types]

sce.dlpfc.tran$cellType_broad <- factor(sce.dlpfc.tran$cellType_broad, levels = exp_cell_types)

table(sce.dlpfc.tran$cellType_broad)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 782        18       388      5455       572      2388      1580

## use ensemblID as rownames
rownames(sce.dlpfc.tran) <- rowData(sce.dlpfc.tran)$gene_id

## save data as sce for downstream compatiabilty 
sce <- sce.dlpfc.tran
save(sce, file = here("processed-data", "12_tran_deconvolution", "sce.dlpfc.tran.Rdata"))

# slurmjobs::job_single('00_data_prep_Tran', create_shell = TRUE, memory = '10G', command = "Rscript 00_data_prep_Tran.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
