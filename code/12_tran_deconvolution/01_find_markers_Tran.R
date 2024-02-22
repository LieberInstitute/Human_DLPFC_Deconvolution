library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("here")
library("sessioninfo")

## load data
here("processed-data", "12_tran_deconvolution", "sce.dlpfc.tran.Rdata")
table(sce$cellType_broad)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     10894      1940     24809     11067

#### Run Marker Finding ####

rowData(sce)$Symbol <- rowData(sce)$Symbol.uniq
rownames(sce) <- rowData(sce)$gene_id

message(Sys.time(), " Find Mean Ratio Genes")
mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_broad", add_symbol = TRUE)

message(Sys.time(), " Find 1vAll genes")
markers_1vALL <- findMarkers_1vAll(sce, cellType_col = "cellType_broad", mod = "~BrNum")

marker_stats <- mean_ratio |> left_join(markers_1vALL)

message(Sys.time(), " Save marker gene data")

data_dir <- here("processed-data", "12_tran_deconvolution", "01_find_markers_Tran")
if(!dir.exists(data_dir)) dir.create(data_dir)
save(marker_stats, file = here("processed-data", "12_tran_deconvolution", "01_find_markers_Tran", "Tran_marker_stats_broad.Rdata"))

# slurmjobs::job_single('01_find_markers_Tran', create_shell = TRUE, memory = '25G', command = "Rscript 01_find_markers_Tran")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
