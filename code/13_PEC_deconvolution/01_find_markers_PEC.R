library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("here")
library("sessioninfo")

## load our DLPFC data
load(here("processed-data", "13_PEC_deconvolution", "sce_PTSDBrainomics.Rdata"), verbose = TRUE)

## drop cell types 
sce <- sce[, sce$cellType_broad != "drop"]
sce$cellType_broad <- droplevels(sce$cellType_broad)

table(sce$cellType_broad)
# Astro EndoMural     Excit     Inhib     Micro     Oligo       OPC 
# 17605       985     86129     36642      9419     35866     10498 

#### Run Marker Finding ####
rowData(sce)$Symbol <- rowData(sce)$gene_name

message(Sys.time(), " Find Mean Ratio Genes")
mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_broad", add_symbol = TRUE)

message(Sys.time(), " Find 1vAll genes")
markers_1vALL <- findMarkers_1vAll(sce, cellType_col = "cellType_broad", mod = "~BrNum")

marker_stats <- mean_ratio |> left_join(markers_1vALL)

message(Sys.time(), " Save marker gene data")
save(marker_stats, file = here("processed-data", "13_PEC_deconvolution", "01_find_markers_PEC", "PEC_marker_stats_broad.Rdata"))

# slurmjobs::job_single('01_find_markers_PEC', create_shell = TRUE, memory = '100G', command = "Rscript 01_find_markers_PEC.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
