
library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
library("here")
library("sessioninfo")

## load our DLPFC data
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

## drop Ambiguous nuc
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

table(sce$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 3979      2157      1601     10894      1940     24809     11067 

# cell_type_colors <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]

#### Run Marker Finding ####

rowData(sce)$Symbol <- rowData(sce)$gene_name
rownames(sce) <- rowData(sce)$gene_id

message(Sys.time(), " Find Mean Ratio Genes")
mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_broad_hc", add_symbol = TRUE)

message(Sys.time(), " Find 1vAll genes")
markers_1vALL <- findMarkers_1vAll(sce, cellType_col = "cellType_broad_hc", mod = "~BrNum")

marker_stats <- mean_ratio |> left_join(markers_1vALL)

message(Sys.time(), " Save marker gene data")
save(marker_stats, file = here("processed-data","06_marker_genes","03_find_markers_broad","marker_stats_broad.Rdata"))

# sgejobs::job_single('03_find_markers_broad', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 03_find_markers_broad.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

