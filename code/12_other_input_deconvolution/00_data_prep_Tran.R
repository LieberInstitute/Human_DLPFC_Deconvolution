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
table(sce.dlpfc.tran$cellType)
# Astro    Excit_A    Excit_B    Excit_C    Excit_D    Excit_E    Excit_F    Inhib_A    Inhib_B    Inhib_C 
# 782        529        773        524        132        187        243        333        454        365 
# Inhib_D    Inhib_E    Inhib_F Macrophage      Micro      Mural      Oligo        OPC      Tcell 
# 413          7          8         10        388         18       5455        572          9

sce.dlpfc.tran$cellType_broad <- gsub("_[A-Z]","",sce.dlpfc.tran$cellType)
sce.dlpfc.tran$cellType_broad <- gsub("Mural", "EndoMural", sce.dlpfc.tran$cellType_broad)
exp_cell_types <- c("Astro", "EndoMural", "Micro", "Oligo","OPC", "Excit", "Inhib") ## cell types considered in this experiment

table(sce.dlpfc.tran$cellType_broad %in% exp_cell_types,
      sce.dlpfc.tran$cellType_broad)
# Astro EndoMural Excit Inhib Macrophage Micro Oligo  OPC Tcell
# FALSE     0         0     0     0         10     0     0    0     9
# TRUE    782        18  2388  1580          0   388  5455  572     0

sce.dlpfc.tran <- sce.dlpfc.tran[,sce.dlpfc.tran$cellType_broad %in% exp_cell_types] ## drop Micro & Tcell
sce.dlpfc.tran$cellType_broad <- factor(sce.dlpfc.tran$cellType_broad, levels = exp_cell_types)
any(is.na(sce.dlpfc.tran$cellType_broad))

table(sce.dlpfc.tran$cellType_broad)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 782        18       388      5455       572      2388      1580

## use ensemblID as rownames
rownames(sce.dlpfc.tran) <- rowData(sce.dlpfc.tran)$gene_id

## save data as sce for downstream compatiabilty 
sce <- sce.dlpfc.tran
save(sce, file = here("processed-data", "12_tran_deconvolution", "sce.dlpfc.tran.Rdata"))

## cell type prop
ct_prop <- colData(sce) |>
  as.data.frame() |>
  dplyr::group_by(donor, cellType_broad) |>
  dplyr::count() |>
  dplyr::group_by(donor) |>
  dplyr::mutate(prop = n/sum(n), 
         dx = "Control")

write.csv(ct_prop, file = here("processed-data", "12_tran_deconvolution", "tran_ct_prop.csv"))

# slurmjobs::job_single('00_data_prep_Tran', create_shell = TRUE, memory = '10G', command = "Rscript 00_data_prep_Tran.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
