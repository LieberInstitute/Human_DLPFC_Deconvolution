library("DeconvoBuddies")
library("SingleCellExperiment")
library("tidyverse")
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

sce <- sce.dlpfc.tran
table(sce$cellType)
colData(sce)

sce$cellType_broad <- gsub("_[A-Z]","",sce$cellType)

ct_tab <- table(sce$cellType_broad)
my_cell_types <- names(ct_tab[ct_tab > 50])
sce <- sce[,sce$cellType_broad %in% my_cell_types]

table(sce.dlpfc.tran$cellType_broad)

## drop Ambiguous nuc
sce <- sce[, sce$cellType_broad != "Ambiguous"]
sce$cellType_broad <- droplevels(sce$cellType_broad)

table(sce$cellType_broad)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     10894      1940     24809     11067

# cell_type_colors <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad)]

#### Run Marker Finding ####

rowData(sce)$Symbol <- rowData(sce)$gene_name
rownames(sce) <- rowData(sce)$gene_id

message(Sys.time(), " Find Mean Ratio Genes")
mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_broad", add_symbol = TRUE)

message(Sys.time(), " Find 1vAll genes")
markers_1vALL <- findMarkers_1vAll(sce, cellType_col = "cellType_broad", mod = "~BrNum")

marker_stats <- mean_ratio |> left_join(markers_1vALL)

message(Sys.time(), " Save marker gene data")
save(marker_stats, file = here("processed-data", "06_marker_genes", "03_find_markers_broad", "marker_stats_broad.Rdata"))

# sgejobs::job_single('03_find_markers_broad', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 03_find_markers_broad.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
